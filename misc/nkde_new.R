# ------------------------------------------------------------
# Network Incident Hotspots via NKDE (spNetwork)
# Author: Natchapon Jongwiriyanurak
# Source pkgs: https://github.com/JeremyGelb/spNetwork
# Data: OSM road network + TRAMS motorcycle-involved incidents
# ------------------------------------------------------------

# ---- packages ----
library(sf)
library(dplyr)
library(spNetwork)
library(units)
library(tmap)
library(future)

# ---- config ----
options(digits = 20)
# setwd(".../NKDE")  # <- set your working directory here

# Parallel notes:
# - On macOS/Linux you can use multicore; on Windows use multisession.
# - spNetwork::nkde.mc and lixelize_lines.mc already leverage parallelism.
# future::plan(multisession, workers = 8)  # optional

# ---- data ----
# Network: must be LINESTRING/MULTILINESTRING in a metric CRS (metres)
network <- st_read("roadnetwork/tha_highway_coverage.gpkg", quiet = TRUE)

# Incidents: POINTS in any CRS; we’ll reproject to match the network
accidents <- st_read("DOH/motorcycle_accidents_TRAMS.shp", quiet = TRUE)
accidents$weight <- accidents$Fatalities*5 + accidents$Serious*3 + accidents$Minor*1

accidents$weight 
# how many accidents$weight =0
sum(accidents$weight ==0)

st_crs(samples_500)

# Optional pre-made midpoints/samples (not required for NKDE itself)
# samples_500 <- st_read("roadnetwork/tha_highway_nkde_samples_500m.gpkg", quiet = TRUE)

# ---- sanity checks / CRS alignment ----
# Ensure network has a projected CRS (metres). If not, pick a suitable one (e.g., UTM 47N: EPSG:32647).
if (is.na(st_crs(network)) || !st_is_longlat(network)) {
  # projected already or CRS missing -> we assume it is projected as per file
  message(sprintf("Network CRS: %s", st_crs(network)$input))
} else {
  stop("Network is in geographic CRS (degrees). Reproject to a metric CRS (e.g., EPSG:32647) before running NKDE.")
}

# Reproject accidents to match network CRS
accidents <- st_transform(accidents, st_crs(network))

# Clean geometries just in case
# network <- st_zm(network, drop = TRUE, what = "ZM")
network <- st_make_valid(network)

# Give network a unique id for later joins
network <- network %>% mutate(net_id = row_number())
network <- st_cast(network)
network$geom <- st_geometry(network)

dens <- nkde.mc(
  lines   = network,
  events  = accidents,   # POINTS; nkde handles snapping internally
  samples = samples_500,     # lixel size in metres (from prepared lixels)
  w       = accidents$weight,        # optional vector of weights (e.g., severity)
  bw      = bw,
  kernel  = "quartic",   # "quartic" or "gaussian" (quartic is common in road safety)
  div     = "lixel",     # density per lixel
  method  = "simple",    # "simple", "continuous", or "discontinuous"
  digits  = 8,
  grid_shape = c(10,10),
  agg = 15,
  sparse = TRUE,
  verbose = TRUE
)

# ---- helper: evaluate HR/PAI at coverage thresholds ----
# coverage_vec is in proportions of total network length (e.g., c(0.05, 0.10, 0.20))
evaluate_hotspots <- function(lixels, density_col = "density",
                              incidents, coverage_vec = c(0.05, 0.10, 0.20),
                              snap_tol = set_units(25, "m")) {
  stopifnot(inherits(lixels, "sf"), inherits(incidents, "sf"))
  # lengths in metres
  lixels <- lixels %>%
    mutate(lx_len_m = as.numeric(st_length(geometry))) %>%
    arrange(desc(.data[[density_col]]))
  
  total_len <- sum(lixels$lx_len_m, na.rm = TRUE)
  if (total_len == 0) stop("Total lixel length is zero after splitting. Check your network.")
  
  # Snap incidents to their nearest lixel (for robust, fast overlay)
  nearest_idx <- st_nearest_feature(incidents, lixels)
  nearest_dist <- st_distance(incidents, lixels[nearest_idx, ], by_element = TRUE)
  
  # Keep only incidents near the network (<= snap_tol)
  keep <- as.numeric(nearest_dist) <= as.numeric(snap_tol)
  inc_keep <- incidents[keep, , drop = FALSE]
  idx_keep <- nearest_idx[keep]
  
  # Total incidents considered for HR
  total_inc <- nrow(inc_keep)
  if (total_inc == 0) {
    warning("No incidents within snapping tolerance of the network; HR = 0 for all coverages.")
  }
  
  # Prepare cumulative length selection by density
  lixels$cum_len <- cumsum(lixels$lx_len_m)
  out_rows <- list()
  
  for (p in coverage_vec) {
    cutoff_len <- p * total_len
    sel <- lixels %>% filter(cum_len <= cutoff_len)
    sel_ids <- sel$lx_id %||% sel$net_id %||% seq_len(nrow(sel)) # fallback if lx_id not present
    
    # Count incidents whose nearest lixel is in selection
    hits <- sum(idx_keep %in% st_drop_geometry(sel) %>% pull(net_id), na.rm = TRUE)
    
    HR <- if (total_inc > 0) hits / total_inc else 0
    PAI <- if (p > 0) HR / p else NA_real_
    
    out_rows[[length(out_rows) + 1]] <- tibble::tibble(
      coverage = p,
      total_incidents = total_inc,
      hits = hits,
      hit_rate = HR,
      pai = PAI,
      sel_len_m = sum(sel$lx_len_m, na.rm = TRUE),
      total_len_m = total_len
    )
  }
  dplyr::bind_rows(out_rows)
}

# ---- sweep across lixel sizes and bandwidths ----
samples_sizes <- c(50, 100, 200, 300, 500)  # metres
bws      <- c(100, 200, 300)           # metres
coverages <- c(0.05, 0.10, 0.20)       # 5%, 10%, 20% of total lixel length

results <- list()

for (samples in samples_sizes) {
  message(sprintf("Lixelising at %d m ...", lx))
  # lixelize_lines.mc splits lines into ~equal segments ("lixels")
  # mindist = lx/2 avoids generating ultra-tiny fragments around nodes
  # lixels <- lixelize_lines.mc(network,
  #                             lx_length = lx,
  #                             mindist = lx / 2,
  #                             verbose = TRUE)
  # 
  # # Give lixels an id; keep parent id if available
  # lixels <- lixels %>%
  #   mutate(net_id = row_number()) %>%
  #   st_make_valid()
  
  for (bw in bws) {
    message(sprintf("  NKDE with bw = %d m ...", bw))
    
    # nkde.mc returns a numeric density vector aligned with 'lixels'
    dens <- nkde.mc(
      lines   = network,
      events  = accidents,   # POINTS; nkde handles snapping internally
      samples = samples,     # lixel size in metres (from prepared lixels)
      w       = rep(1,nrow(accidents$weight)),        # optional vector of weights (e.g., severity)
      bw      = bw,
      kernel  = "quartic",   # "quartic" or "gaussian" (quartic is common in road safety)
      div     = "lixel",     # density per lixel
      method  = "simple",    # "simple", "continuous", or "discontinuous"
      digits  = 8,
      verbose = TRUE
    )
    
    # Attach density to lixels
    lixels$density <- as.numeric(dens)
    
    # Evaluate HR/PAI at requested coverages
    eval_tbl <- evaluate_hotspots(
      lixels   = lixels %>% mutate(lx_id = row_number()),
      density_col = "density",
      incidents   = accidents,
      coverage_vec = coverages,
      snap_tol     = set_units(25, "m")
    ) %>%
      mutate(lixel_size_m = samples,
             bandwidth_m  = bw)
    
    results[[paste0("lx", lx, "_bw", bw)]] <- eval_tbl
    
    # Optional: save hotspot layer for mapping
    hotspot_out <- lixels %>% arrange(desc(density))
    st_write(hotspot_out,
             sprintf("out/hotspots_lx%03d_bw%03d.gpkg", samples, bw),
             layer = "hotspots",
             delete_layer = TRUE)
  }
}

# ---- combine & view results ----
results_tbl <- dplyr::bind_rows(results) %>%
  arrange(lixel_size_m, bandwidth_m, coverage)

print(results_tbl)

# Save a tidy CSV for the paper
# write.csv(results_tbl, "out/hotspot_eval_hr_pai.csv", row.names = FALSE)

# ---- quick map (optional): top 10% by length for a given combo ----
# Choose one combo to visualise
viz_lx <- 100
viz_bw <- 200

# Recompute density for that combo (or reuse from loop if kept)
lixels_v <- lixelize_lines.mc(network, lx_length = viz_lx, mindist = viz_lx/2)
lixels_v$density <- nkde.mc(
  lines = lixels_v, events = accidents,
  bw = viz_bw, kernel = "quartic", div = "lixel", method = "simple", verbose = TRUE
)
lixels_v <- lixels_v %>% mutate(lx_len_m = as.numeric(st_length(geometry))) %>%
  arrange(desc(density))
tot_len <- sum(lixels_v$lx_len_m)
sel10 <- lixels_v %>% mutate(cum_len = cumsum(lx_len_m)) %>% filter(cum_len <= 0.10 * tot_len)

tmap_mode("view")
tm_shape(lixels_v) + tm_lines(col = "grey80") +
  tm_shape(lixels_v) + tm_lines(col = "density", palette = "magma", lwd = 3, legend.col.show = TRUE) +
  tm_shape(sel10)    + tm_lines(col = "turquoise", lwd = 6) +
  tm_shape(accidents) + tm_dots(col = "red", size = 0.05) +
  tm_view(set.zoom.limits = c(7, 17))
