# This Script used to identify Network incident Hotspots
# The process spNetwork package (https://github.com/JeremyGelb/spNetwork) by Gelb Jérémy (2021)
# The data used in this script includes: OSM network and Motorcycle-involved data from TRAMS

# Load required libraries
library(spNetwork)
library(tmap)
library(sf)
library(dplyr)
library(future)


# Set the working directory and locale
setwd("/Users/stupong/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Natchapon PhD/Upgrading/PhD/NKDE")
# Sys.setlocale(locale = "Thai")
options(digits = 20)

# nbrOfWorkers()
# nbrOfFreeWorkers()
# # Set the multisession plan
# future::plan(future::multisession(workers = 8))
# if (!inherits(future::plan(), "sequential")) {
#   future::plan(future::sequential)
# }

# Read network and transform CRS


network <- st_read("roadnetwork/thailand-250902-free.shp/gis_osm_roads_free_1.shp")
network_highway <- st_read("roadnetwork/tha_highway_coverage_clip.gpkg", layer="coverage")

network <- st_make_valid(network) |> st_zm(drop = TRUE, what = "ZM")
network <- st_collection_extract(network, "LINESTRING")
network <- st_cast(network, "LINESTRING", warn = FALSE)

network$len_m <- as.numeric(st_length(network))
unique(st_geometry_type(network))  # should be only "LINESTRING"

summary(as.numeric(st_length(network)))
hist(as.numeric(st_length(network)), breaks = 50)


lixels_100 <- lixelize_lines.mc(
  network,
  lx_length = 100,
  mindist   = 50,
  verbose   = TRUE
)


lixels_100 <- lixelize_lines.mc(
  net,
  lx_length = 100,   # 100 m lixels
  mindist   = 50,    # >= 50 m between sample points on the *same* line
  verbose   = TRUE
)


plot(network[1:1000,])

samples <- lines_points_along(network, 100)

# lixels_50 <- lixelize_lines(network, lx_length  = 50, verbose=TRUE)
# samepls_50 <- sample_lines.mc(lixels, n = 1)

lixels_100 <- lixelize_lines.mc(network, lx_length = 100, mindist = 50, verbose=TRUE)
plot(lixels_100[1:1000,])
samepls_100 <- sample_lines.mc(lixels,verbose=TRUE)

lixels_200 <- lixelize_lines.mc(network, lx_length = 200)
samepls_200 <- sample_lines.mc(lixels, n = 1, progress = TRUE)

lixels_300 <- lixelize_lines.mc(network, lx_length = 300, verbose=TRUE)
samepls_300 <- sample_lines.mc(lixels, n = 1, progress = TRUE)


#save lixels and samples as shp
st_write(lixels_100, "lixels_100.gpkg")
st_write(samepls_100, "samepls_100.gpkg")
st_write(lixels_200, "lixels_200.gpkg")
st_write(samples_200, "sameples_200.gpkg")
st_write(lixels_300, "lixels_300.gpkg")
st_write(samples_300, "sameples_300.gpkg")

# Read TRAMS
accidents <- st_read("DOH/motorcycle_accidents_TRAMS.shp")

# plot [:1000] of network and accidents
tm_shape(network) + tm_lines() +
  tm_shape(accidents[1:1000,]) + tm_dots(col = "red", size = 0.1)




# 1) Build lixels for each size
lx_sizes <- c(50, 100, 200, 300, 500)
bws <- c(100, 200, 300)
out <- list()

for (lx in lx_sizes) {
  lixels <- lixelize_lines.mc(network, lx_length = lx, mindist = lx/2)
  # 2) NKDE for each bandwidth
  for (bw in bws) {
    kd <- nkde.mc(
      lines = lixels,
      events = crashes,               # sf POINTS snapped or raw (nkde handles snapping)
      w = NULL,                       # or severity weights
      bw = bw,
      kernel = "quartic",             # or "gaussian"
      div = "lixel",                  # density per lixel
      method = "simple"               # or "discontinuous"/"continuous" as available
    )
    # 3) Rank lixels by density, compute HR/PAI at 5/10/20%
    res <- evaluate_hotspots(kd, crashes, lixels, coverages = c(0.01, 0.05, 0.10))
    out[[paste0("lx",lx,"_bw",bw)]] <- res
  }
}

# evaluate_hotspots should:
# - sort lixels by density desc
# - for each coverage p: select top p% length, count crashes within -> HR = Cp/CT
# - PAI = HR / p


