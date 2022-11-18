
library(raster)
library(data.table)
library(posim)

setwd("/Volumes/xiaopenyou/posim")
proj4 <- "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs"

rasters = list(
    bathy = raster("bathy.asc", crs = proj4),
    blocks = raster("blocks.asc", crs = proj4),
    abundance = raster("abundance.asc", crs = proj4),
    patches = raster("patches.asc", crs = proj4),
    disttocoast = raster("disttocoast.asc", crs = proj4),
    fishery = raster("fishery.asc", crs = proj4),
    quarter1 = raster("quarter1.asc", crs = proj4),
    quarter2 = raster("quarter2.asc", crs = proj4),
    quarter3 = raster("quarter3.asc", crs = proj4),
    quarter4 = raster("quarter4.asc", crs = proj4)
)

rasters <- lapply(rasters, crop, extent(350000, 600000, 7400000, 7750000))

# recalibrate

cells <- which(rasters$bathy[] < 0)
patch_pool <- sample(cells, size = length(cells))

patch_count <- 9200
rasters$patches[] <- 0
rasters$patches[patch_pool[1:patch_count]] <- 1

sasc(filename = "/Volumes/xiaopenyou/posim/Lofoten.sasc",
     bathymetry = rasters$bathy,
     distToCoast = rasters$disttocoast,
     block = rasters$blocks,
     abundanceRegion = rasters$abundance,
     fishery = rasters$fishery,
     maxent1 = rasters$quarter1,
     maxent2 = rasters$quarter2,
     maxent3 = rasters$quarter3,
     maxent4 = rasters$quarter4,
     foodLevel = rasters$patches)


x <- posim(steps = 3*365*48,
         follow = 10,
         start = 1, 
         sasc = "/Volumes/xiaopenyou/posim/Lofoten.sasc",
         distEnergyMultiplier = 0,
         N = 100)
plot(x)
plot(x, track = TRUE, select = 1)
plot(x, track = TRUE, highlight = 1)
animate(x, steps = 365*48, filename = "tracks.gif")


library(sf)
ext <- extent(rasters$bathy)
xstep <- (ext[2]-ext[1]) / ncol(rasters$bathy)
ystep <- (ext[4]-ext[3]) / nrow(rasters$bathy)
track <- as.data.table(x@tracks)
track[, x := ext[1] + x * xstep]
track[, y := ext[3] + y * ystep]
track <- st_as_sf(track, coords = c("x", "y"), crs = 32633)
plot(rasters$bathy)
plot(head(track[track$id == 3,],n  =48*365), add=TRUE, type="l")

xy <- cbind(xyF(st_coordinates(track[track$id == 3,]))

library(ggplot2)
library(rasterVis)
gplot(rasters$bathy, maxpixels = ncell(rasters$bathy)) +
    geom_tile(aes(fill=value)) +
    geom_path(data = xy, aes(x = X, y = Y))+
    coord_equal(expand=FALSE)


landscape <- cbind(xyFromCell(rasters$bathy, 1:ncell(rasters$bathy)), getValues(rasters$bathy))
landscape <- as.data.frame(landscape)
names(landscape) <- c("x", "y", "Bathymetry")

ggplot(data = landscape) +
    geom_tile(aes(x, y, fill = Bathymetry)) +
    coord_equal(expand = FALSE)

