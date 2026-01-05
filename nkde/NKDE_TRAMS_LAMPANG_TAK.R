# ------------------------------------------------------------
# NKDE for two provinces (LAMPANG, TAK) with tile auto-repair
# Keeps ALL columns from sample layers and adds `density`
# Requires: sf, dplyr, spNetwork, stringr
# ------------------------------------------------------------

library(sf)
library(dplyr)
library(spNetwork)
library(stringr)

# ---------- INPUTS ----------
# If not already in memory, read them (edit paths as needed):
# network   <- st_read("roadnetwork/tha_highway_coverage.gpkg", quiet = TRUE)
# accidents <- st_read("DOH/motorcycle_accidents_TRAMS.shp", quiet = TRUE)
# provinces <- st_read(".../TH_Province.shp", quiet = TRUE)
# samples_50  <- st_read("roadnetwork/tha_highway_nkde_samples_50m.gpkg",  quiet=TRUE)
# samples_100 <- st_read("roadnetwork/tha_highway_nkde_samples_100m.gpkg", quiet=TRUE)
# samples_200 <- st_read("roadnetwork/tha_highway_nkde_samples_200m.gpkg", quiet=TRUE)
# samples_300 <- st_read("roadnetwork/tha_highway_nkde_samples_300m.gpkg", quiet=TRUE)
# samples_500 <- st_read("roadnetwork/tha_highway_nkde_samples_500m.gpkg", quiet=TRUE)

# Map spacings to objects for iteration:
sample_map <- list(
  `50`  = samples_50,
  `100` = samples_100,
  `200` = samples_200,
  `300` = samples_300,
  `500` = samples_500
)

# Bandwidths and spacings to run:
bws         <- c(100, 200, 300)          # metres
spacing_vec <- c(50, 100, 200, 300, 500) # metres

# Output directory:
out_dir <- "out/nkde_lampang_tak_new"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- BASIC CHECKS & CRS ----------
stopifnot(!st_is_longlat(network))                 # must be projected (metres)
accidents <- st_transform(accidents, st_crs(network))
provinces <- st_transform(provinces, st_crs(network))
sample_map <- lapply(sample_map, \(s) st_transform(s, st_crs(network)))

# Robust accident weights (fallback = 1)
if (!("weight" %in% names(accidents))) {
  f <- if ("Fatalities" %in% names(accidents)) accidents$Fatalities else 0
  s <- if ("SeriousInj"  %in% names(accidents)) accidents$SeriousInj  else 0
  m <- if ("MinorInjur"  %in% names(accidents)) accidents$MinorInjur  else 0
  w <- (dplyr::coalesce(f,0) * 6) + (dplyr::coalesce(s,0) * 4) + (dplyr::coalesce(m,0) * 2)
  w[!is.finite(w) | is.na(w) | w <= 0] <- 1
  accidents$weight <- w
}

# ---------- HELPERS ----------
# simple grid over polygon (nx * ny)
make_tiles <- function(poly_sf, nx = 3, ny = 3) {
  bb  <- st_bbox(poly_sf)
  xs  <- seq(bb["xmin"], bb["xmax"], length.out = nx + 1)
  ys  <- seq(bb["ymin"], bb["ymax"], length.out = ny + 1)
  polys <- vector("list", nx * ny); k <- 1
  for (i in seq_len(nx)) for (j in seq_len(ny)) {
    polys[[k]] <- st_polygon(list(rbind(
      c(xs[i], ys[j]), c(xs[i+1], ys[j]),
      c(xs[i+1], ys[j+1]), c(xs[i], ys[j+1]),
      c(xs[i], ys[j])
    )))
    k <- k + 1
  }
  grid <- st_sf(geometry = st_sfc(polys, crs = st_crs(poly_sf)))
  st_intersection(grid, st_geometry(poly_sf)) |>
    st_make_valid()
}

# Run NKDE in a single tile; attempt auto-repair if it fails.
# NOTE: We KEEP all columns from the sample subset `smp_t` and only add `density`.
run_tile <- function(tile, net_sub, acc_sub, smp_buf, bw = 100) {
  halo <- 100
  th <- st_buffer(tile, halo)
  
  # Intersections keep attributes; we only clean geometry type where needed
  net_t <- st_intersection(net_sub, th) |>
    st_collection_extract("LINESTRING", warn = FALSE)
  acc_t <- st_intersection(acc_sub, th) |>
    st_collection_extract("POINT", warn = FALSE)
  smp_t <- st_intersection(smp_buf, th)  |>
    st_collection_extract("POINT", warn = FALSE)
  
  if (nrow(smp_t) == 0 || nrow(net_t) == 0) return(list(ok = TRUE, smp = smp_t)) # no samples or no network
  if (nrow(acc_t) == 0) { smp_t$density <- 0; return(list(ok = TRUE, smp = smp_t)) }
  
  wloc <- if ("weight" %in% names(acc_t)) acc_t$weight else rep(1, nrow(acc_t))
  
  out <- try(nkde.mc(
    lines   = net_t,
    events  = acc_t,
    samples = smp_t,
    w       = wloc,
    bw      = bw,
    kernel  = "quartic",
    div     = "bw",        # IMPORTANT for point samples
    method  = "simple",
    digits  = NA,          # avoid rounding/coercion issues
    sparse  = FALSE,       # return full vector
    verbose = FALSE
  ), silent = TRUE)
  
  if (!inherits(out, "try-error")) {
    smp_t$density <- as.numeric(out)
    return(list(ok = TRUE, smp = smp_t))
  }
  
  # Auto-repair: try dropping a few edges in this tile if it fails
  n_edges <- nrow(net_t)
  cap <- min(n_edges, 60L)
  for (j in seq_len(cap)) {
    out2 <- try(nkde.mc(
      lines   = net_t[-j, ],
      events  = acc_t,
      samples = smp_t,
      w       = wloc,
      bw      = bw,
      kernel  = "quartic",
      div     = "bw",
      method  = "simple",
      digits  = NA,
      sparse  = FALSE,
      verbose = FALSE
    ), silent = TRUE)
    if (!inherits(out2, "try-error")) {
      message(sprintf("      repaired: dropped edge %d/%d", j, n_edges))
      smp_t$density <- as.numeric(out2)
      return(list(ok = TRUE, smp = smp_t))
    }
  }
  
  warning("      tile could not be repaired; assigning zeros")
  smp_t$density <- 0
  list(ok = FALSE, smp = smp_t)
}

# Run a whole province (by name) across spacings and bandwidths
run_province <- function(prov_name,
                         spacing_vec,
                         bws,
                         sample_map,
                         out_dir,
                         grid_nx = 3, grid_ny = 3,
                         prov_buffer = 600) {
  
  prov <- provinces[provinces$PROV_NAME == prov_name, , drop = FALSE]
  stopifnot(nrow(prov) == 1)
  message(sprintf("=== Province: %s ===", prov_name))
  
  prov_buf    <- st_buffer(st_geometry(prov), dist = prov_buffer)
  prov_buf_sf <- st_sf(geometry = prov_buf, crs = st_crs(prov))
  
  # subset inputs to buffered province (attributes preserved)
  net_sub <- st_filter(network,   prov_buf_sf) |>
    st_make_valid() |>
    st_collection_extract("LINESTRING", warn = FALSE)
  acc_sub <- st_filter(accidents, prov_buf_sf) |>
    st_make_valid() |>
    st_collection_extract("POINT", warn = FALSE)
  
  # tiles over the UNBUFFERED province
  tile_sf <- make_tiles(prov, nx = grid_nx, ny = grid_ny)
  
  for (sx in spacing_vec) {
    smp <- sample_map[[as.character(sx)]]
    if (is.null(smp)) {
      message(sprintf("  spacing %d m: samples not found; skipping", sx))
      next
    }
    smp <- st_transform(smp, st_crs(prov)) # ensure CRS match
    
    # samples within buffered province; keep their original attributes
    smp_buf <- st_filter(smp, prov_buf_sf)
    
    for (bw in bws) {
      message(sprintf("  -> spacing = %dm, bw = %dm", sx, bw))
      res_list <- vector("list", length = nrow(tile_sf))
      ok_flag  <- logical(nrow(tile_sf))
      
      for (i in seq_len(nrow(tile_sf))) {
        res <- run_tile(tile_sf[i, ], net_sub, acc_sub, smp_buf, bw = bw)
        ok_flag[i] <- res$ok
        res_list[[i]] <- res$smp  # contains all original sample cols + density (added inside run_tile)
      }
      
      # stitch + keep only samples INSIDE the unbuffered province
      smp_all <- do.call(rbind, res_list)
      if (!is.null(smp_all) && nrow(smp_all)) {
        smp_out <- smp_all[st_intersects(smp_all, prov, sparse = FALSE), , drop = FALSE]
        safe_name <- str_squish(gsub("[^[:alnum:]_]+", "_", prov_name))
        out_file  <- file.path(out_dir, sprintf("nkde_%s_lx%03d_bw%03d.gpkg", safe_name, sx, bw))
        st_write(smp_out, out_file, layer = "nkde_density", delete_layer = TRUE, quiet = TRUE)
        message(sprintf("     saved: %s", out_file))
      } else {
        message("     no samples in province after tiling; nothing saved")
      }
      
      bad_after <- which(!ok_flag)
      if (length(bad_after)) {
        message(sprintf("     tiles still problematic after auto-repair: %s",
                        paste(bad_after, collapse = ", ")))
      }
    }
  }
}

# ---------- RUN: Lampang & Tak ----------
run_province("LAMPANG", spacing_vec, bws, sample_map, out_dir)
run_province("TAK",      spacing_vec, bws, sample_map, out_dir)
