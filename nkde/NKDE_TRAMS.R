# ------------------------------------------------------------
# Network Incident Hotspots via NKDE (spNetwork)
# Author: Natchapon Jongwiriyanurak
# Purpose: Compute NKDE densities for multiple bandwidths (bw)
#          and pre-made sample spacings (midpoints), saving per-province outputs.
# Notes:
#   - We DO NOT create lixels here (you already have samples).
#   - We DO NOT compute HR/PAI here (do later).
#   - We reduce compute by processing province-by-province with a buffer
#     (edge effect handling: include nearby events/network, but only save
#      densities for samples inside the unbuffered province).
# ------------------------------------------------------------

# ---- packages ----
library(sf)
library(dplyr)
library(spNetwork)
library(units)
library(tidyr)
library(stringr)

# ---- config ----
options(digits = 20)
setwd("/Users/stupong/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Natchapon PhD/Upgrading/PhD/NKDE")

# Output folder
out_dir <- "out/nkde"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Bandwidths (metres) and the sample spacings (metres) you already prepared
bws          <- c(100, 200, 300)
sample_sizes <- c(50, 100, 200, 300, 500)

# A safe NKDE config (use dense output to avoid zero-length tile joins)
nkde_kernel  <- "quartic"
nkde_method  <- "simple"
nkde_sparse  <- FALSE         # keep zeros so joins don't fail
nkde_digits  <- 8

# Province buffer (metres) for edge handling:
# Use at least 2 * max(bw) so kernels can "see" across borders.
prov_buffer_m <- 2 * max(bws)

# ---- data paths ----
# Network (metric CRS required)
network_path   <- "roadnetwork/tha_highway_coverage.gpkg"

# Incidents (POINTS). Must contain (or we will make) a 'weight' column
accidents_path <- "DOH/motorcycle_accidents_TRAMS.shp"

# Province boundaries (POLYGON/MULTIPOLYGON)
provinces_path <- "/Users/stupong/Library/CloudStorage/OneDrive-SharedLibraries-UniversityCollegeLondon(2)/SpaceTimeLab - Natchapon PhD - Natchapon PhD/Upgrading/PhD/traffic_incident_analysis/data/boundary/TH_Province.shp"

# Sample files (POINTS along the network at fixed spacing).
# Adapt these filenames if yours differ.
samples_files <- list(
  `50`  = "roadnetwork/tha_highway_nkde_samples_50m.gpkg",
  `100` = "roadnetwork/tha_highway_nkde_samples_100m.gpkg",
  `200` = "roadnetwork/tha_highway_nkde_samples_200m.gpkg",
  `300` = "roadnetwork/tha_highway_nkde_samples_300m.gpkg",
  `500` = "roadnetwork/tha_highway_nkde_samples_500m.gpkg"
)

# ---- read data ----
message("Reading network / accidents / provinces / samples ...")

network   <- st_read(network_path, quiet = TRUE)
accidents <- st_read(accidents_path, quiet = TRUE)
provs     <- st_read(provinces_path, quiet = TRUE)

# Read sample sets into a named list keyed by spacing
samples_list <- lapply(samples_files, \(p) st_read(p, quiet = TRUE))

# ---- CRS alignment & basic cleaning ----
# Ensure everything is in a metric CRS (metres). If the network layer is in UTM, use that as the target.
if (is.na(st_crs(network)) || st_is_longlat(network)) {
  stop("Network must be in a projected (metric) CRS. Reproject before proceeding (e.g., EPSG:32647).")
}

# Reproject others to network CRS
accidents <- st_transform(accidents, st_crs(network))
provs     <- st_transform(provs, st_crs(network))
samples_list <- lapply(samples_list, \(s) st_transform(s, st_crs(network)))

# Clean geometries
network <- network |>
  st_make_valid() |>
  st_cast("LINESTRING", warn = FALSE) |>
  filter(!st_is_empty(geometry))

accidents <- accidents |>
  st_make_valid() |>
  filter(!st_is_empty(geometry))

provs <- provs |>
  st_make_valid() |>
  filter(!st_is_empty(geometry))

samples_list <- lapply(samples_list, \(s) {
  s |>
    st_make_valid() |>
    filter(!st_is_empty(geometry))
})

# ---- weights for accidents ----
# Create a robust weight vector:
# If you have severity fields (Fatalities/Serious/Minor), use them; else weight=1.
accidents$weight <- accidents$Fatalities * 6 +
    accidents$SeriousInj    * 4 +
    accidents$MinorInjur      * 2

#if accidents$weight is 0 or NA, set to 1
accidents$weight[!is.finite(accidents$weight) | accidents$weight <= 0] <- 1

max(accidents$weight)
min(accidents$weight)

# ---- helper: run NKDE for one province, one sample set, one bw ----
run_nkde_one <- function(prov_row, samples_sf, bw) {
  # prov_row: single-row sf polygon of the province (unbuffered)
  # samples_sf: sample points (full country) for a given spacing
  # bw: numeric bandwidth (metres)
  # Returns: sf with samples inside the province and a 'density' column
  
  # Buffered polygon for input selection (edge handling)
  prov_buf <- st_buffer(st_geometry(prov_row), dist = prov_buffer_m)
  prov_buf <- st_as_sf(prov_buf) |> st_set_crs(st_crs(prov_row))
  
  # Clip/select inputs to the buffered province
  net_sub <- network[st_intersects(network, prov_buf, sparse = FALSE), , drop = FALSE]
  if (nrow(net_sub) == 0) return(NULL)
  
  acc_sub <- accidents[st_intersects(accidents, prov_buf, sparse = FALSE), , drop = FALSE]
  if (nrow(acc_sub) == 0) {
    # No events in buffered area: return zero density for samples inside province
    smp_in_prov <- samples_sf[st_intersects(samples_sf, prov_row, sparse = FALSE), , drop = FALSE]
    if (nrow(smp_in_prov) == 0) return(NULL)
    smp_in_prov$density <- 0
    return(smp_in_prov)
  }
  
  # Samples to evaluate density AT (we will compute at samples in buffered area,
  # then keep only samples that fall inside the unbuffered province for output)
  smp_buf  <- samples_sf[st_intersects(samples_sf, prov_buf, sparse = FALSE), , drop = FALSE]
  if (nrow(smp_buf) == 0) return(NULL)
  
  # Weights aligned to acc_sub
  w <- acc_sub$weight
  if (is.null(w)) w <- rep(1, nrow(acc_sub))
  w[!is.finite(w) | is.na(w)] <- 1
  
  # Compute NKDE (safe settings: dense output, no tiling while you debug)
  dens <- nkde.mc(
    lines    = net_sub,      # network subset (LINESTRING)
    events   = acc_sub,      # POINTS within buffered province
    samples  = smp_buf,      # POINTS (midpoints) within buffered province
    w        = w,            # weights for events
    bw       = bw,
    kernel   = nkde_kernel,  # "quartic"
    div      = "lixel",      # per-segment density proxy; ok to use with point samples
    method   = nkde_method,  # "simple"
    digits   = nkde_digits,
    sparse   = nkde_sparse,  # keep zeroes to avoid join-length errors
    verbose  = TRUE
  )
  
  smp_buf$density <- as.numeric(dens)
  
  # Keep only samples INSIDE the (unbuffered) province for saving
  smp_out <- smp_buf[st_intersects(smp_buf, prov_row, sparse = FALSE), , drop = FALSE]
  if (nrow(smp_out) == 0) return(NULL)
  
  smp_out
}

# ---- main loop: per-province, per-sample-size, per-bw ----
# You can subset provinces here if you want to test on one or two first.
# provs <- provs %>% slice(1:2)

name_col <- intersect(names(provs), c("NAME_1", "prov_name", "PROV_NAME", "Name", "ADM1_TH", "ADM1_EN"))
if (length(name_col) == 0) name_col <- names(provs)[1]  # fallback to first column
name_col <- name_col[1]

for (pi in seq_len(nrow(provs))) {
  prov <- provs[pi, , drop = FALSE]
  prov_name <- as.character(prov[[name_col]])
  safe_name <- str_squish(gsub("[^[:alnum:]_]+", "_", prov_name))
  message(sprintf("=== Province %d/%d: %s ===", pi, nrow(provs), prov_name))
  
  for (sx in sample_sizes) {
    samples_sf <- samples_list[[as.character(sx)]]
    if (is.null(samples_sf)) {
      warning(sprintf("No samples for spacing %d m; skipping.", sx))
      next
    }
    
    for (bw in bws) {
      message(sprintf("  -> spacing = %dm, bw = %dm", sx, bw))
      
      smp_den <- run_nkde_one(prov_row = prov, samples_sf = samples_sf, bw = bw)
      if (is.null(smp_den)) {
        message("     (no data in this slice; skipping)")
        next
      }
      
      # Save per-province-per-combo as GPKG
      out_file  <- file.path(out_dir, sprintf("nkde_%s_lx%03d_bw%03d.gpkg", safe_name, sx, bw))
      out_layer <- "nkde_density"
      
      # Overwrite the layer if it exists (file may be reused across combos; we want one layer per combo)
      st_write(smp_den, out_file, layer = out_layer, delete_layer = TRUE, quiet = TRUE)
      message(sprintf("     saved: %s", out_file))
    }
  }
}

message("Done. Densities saved per province & parameter combo.")

# ---- helper: run NKDE for one province, one sample set, one bw ----
run_nkde_one <- function(prov_row, samples_sf, bw) {
  # Buffered polygon for input selection (edge handling)
  prov_buf <- st_buffer(st_geometry(prov_row), dist = prov_buffer_m)
  prov_buf <- st_as_sf(prov_buf) |> st_set_crs(st_crs(prov_row))
  
  # Subset inputs to buffer
  net_sub <- network[st_intersects(network, prov_buf, sparse = FALSE), , drop = FALSE]
  if (nrow(net_sub) == 0) return(NULL)
  
  acc_sub <- accidents[st_intersects(accidents, prov_buf, sparse = FALSE), , drop = FALSE]
  # If no events in buffer, return zeros for samples inside province
  smp_in_prov <- samples_sf[st_intersects(samples_sf, prov_row, sparse = FALSE), , drop = FALSE]
  if (nrow(acc_sub) == 0) {
    if (nrow(smp_in_prov) == 0) return(NULL)
    smp_in_prov$density <- 0
    return(smp_in_prov)
  }
  
  # Samples used for computation (buffered), then we’ll keep only those inside the province
  smp_buf <- samples_sf[st_intersects(samples_sf, prov_buf, sparse = FALSE), , drop = FALSE]
  if (nrow(smp_buf) == 0) return(NULL)
  
  # Robust weights
  w <- acc_sub$weight
  if (is.null(w)) w <- rep(1, nrow(acc_sub))
  w[!is.finite(w) | is.na(w) | w <= 0] <- 1
  
  # ---- NKDE: IMPORTANT: div = "bw" when samples are POINTS ----
  dens <- nkde.mc(
    lines    = net_sub,
    events   = acc_sub,
    samples  = smp_buf,     # POINTS
    w        = w,
    bw       = bw,
    kernel   = nkde_kernel, # "quartic"
    div      = "bw",        # <— change here
    method   = nkde_method, # "simple" is fine
    digits   = nkde_digits,
    sparse   = nkde_sparse, # keep zeros so lengths match
    verbose  = TRUE
  )
  
  smp_buf$density <- as.numeric(dens)
  
  # Keep only samples INSIDE the (unbuffered) province for output
  smp_out <- smp_buf[st_intersects(smp_buf, prov_row, sparse = FALSE), , drop = FALSE]
  if (nrow(smp_out) == 0) return(NULL)
  
  smp_out
}


# ---- resume from province 41 ---- error (Lampang 41, TAK 51) 38 บึงกาฬ NONG KHAI
provs <- provinces[27, ]   # uncomment to start from 41
provs$PROV_NAME <- "Bueng Kan"

provs


# run_nkde_one <- function(prov_row, samples_sf, bw) {
#   # Buffered polygon for input selection (edge handling)
#   prov_buf <- st_buffer(st_geometry(prov_row), dist = prov_buffer_m)
#   prov_buf <- st_as_sf(prov_buf) |> st_set_crs(st_crs(prov_row))
#   
#   # Subset network to buffer
#   net_sub <- network[st_intersects(network, prov_buf, sparse = FALSE), , drop = FALSE]
#   if (nrow(net_sub) == 0) return(NULL)
#   
#   # ---- SANITISE LINES ----
#   net_sub <- net_sub |>
#     st_make_valid() |>
#     st_collection_extract("LINESTRING", warn = FALSE)
#   # drop empties, zero-length, or 1-node lines
#   net_sub <- net_sub[!st_is_empty(net_sub), , drop = FALSE]
#   net_sub <- net_sub[as.numeric(st_length(net_sub)) > 0, , drop = FALSE]
#   net_sub <- net_sub[st_npoints(net_sub) >= 2, , drop = FALSE]
#   if (nrow(net_sub) == 0) return(NULL)
#   
#   # Avoid name clash with spNetwork’s default edge weight "length"
#   if ("length" %in% names(net_sub)) {
#     names(net_sub)[names(net_sub) == "length"] <- "length_attr"
#   }
#   
#   # Subset incidents & samples to buffer
#   acc_sub <- accidents[st_intersects(accidents, prov_buf, sparse = FALSE), , drop = FALSE]
#   smp_in_prov <- samples_sf[st_intersects(samples_sf, prov_row,  sparse = FALSE), , drop = FALSE]
#   if (nrow(acc_sub) == 0) {
#     if (nrow(smp_in_prov) == 0) return(NULL)
#     smp_in_prov$density <- 0
#     return(smp_in_prov)
#   }
#   smp_buf <- samples_sf[st_intersects(samples_sf, prov_buf, sparse = FALSE), , drop = FALSE]
#   if (nrow(smp_buf) == 0) return(NULL)
#   
#   # Robust weights
#   w <- acc_sub$weight
#   if (is.null(w)) w <- rep(1, nrow(acc_sub))
#   w[!is.finite(w) | is.na(w) | w <= 0] <- 1
#   
#   # NKDE at POINT samples -> use div="bw"
#   dens <- nkde.mc(
#     lines    = net_sub,
#     events   = acc_sub,
#     samples  = smp_buf,
#     w        = w,
#     bw       = bw,
#     kernel   = nkde_kernel,  # "quartic"
#     div      = "bw",         # IMPORTANT with POINT samples
#     method   = nkde_method,  # "simple"
#     digits   = 6,            # slightly less rounding to avoid vertex collapse
#     sparse   = nkde_sparse,  # keep zeros
#     verbose  = TRUE
#   )
#   
#   smp_buf$density <- as.numeric(dens)
#   
#   # Keep only samples INSIDE the unbuffered province
#   smp_out <- smp_buf[st_intersects(smp_buf, prov_row, sparse = FALSE), , drop = FALSE]
#   if (nrow(smp_out) == 0) return(NULL)
#   
#   smp_out
# }






