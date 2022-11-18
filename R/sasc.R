#' Creates a serialized (flattened) raster stack from raster files
#' 
#' Creates a sasc data file from raster files, optionally cropping raster data, that can be used by [posim::sim()]
#' The so-called serialized ASCII (sasc) file is a sort of flattened raster stack, where the cell values
#' of different rasters are used as column vectors in a data file. The first line in the sasc file is a metadata header
#' specifying the extent of the landscape, block size, number of food patches and the max values for each of
#' the different columns (for more convenient processing when reading data again later). 
#' Note that all rasters used must have the same extent and resolution.
#' 
#' @rdname sasc
#' @param filename Output filename
#' @param bathymetry Raster. Average depths.
#' @param distToCoast Raster. Distance (in units of cells) to the nearest point on the coast
#' @param abundanceRegion Raster. Abundance region (if missing, the whole landscape is grouped into one abundance region)
#' @param fisheryRegion Raster. fishery region (if missing, the whole landscape is grouped into one fishery region)
#' @param block Raster. Block membership of each cell. Alternatively, use block_size to automatically assign blocks
#' @param foodLevel Raster. Denotes whether a cell is a food patch or not (1 = food patch, 0 otherwise)
#' @param block_size Single integer. The size of blocks in units of cells. Set this to automatically divide the landscape into blocks. Ignored if \emph{block} is specified
#' @returns None
#' @md
#' @export
sasc <- function(filename, bathymetry, distToCoast, block = NULL, block_size = NULL, abundanceRegion = NULL, fisheryRegion = NULL, foodLevel, 
                 maxent1, maxent2 = NULL, maxent3 = NULL, maxent4 = NULL, precision = 2) {
    
    scipen = options("scipen")
    options(scipen = 100000)
    
    if (is.null(block)) {
        if (is.null(block_size)) {
            stop("You must provide either block or block_size")
        } else if (block_size > ncol(r) | block_size > nrow(r))  {
            stop("block_size cannot be larger than the number of rows/columns")
        } else {
            block <- posim::blockify(bathymetry, block_size)
        }
    }
    
    block <- as.integer(raster::values(block))
    distToCoast <- raster::values(distToCoast)
    foodLevel <- raster::values(foodLevel)
    
    ncol <- ncol(bathymetry)
    nrow <- nrow(bathymetry)
    
    bathymetry <- raster::values(bathymetry)
    
    if (all(bathymetry >= 0, na.rm = TRUE)) { # if there are no negative vals, then depth must have been specified as positive values
        bathymetry <- -1 * bathymetry
        warning("Bathymetry automatically converted to negative values");
    }
    
    water_cells <- which(!is.na(bathymetry) & bathymetry < 0)
    
    if (is.null(maxent2)) maxent2 = maxent1
    if (is.null(maxent3)) maxent3 = maxent1
    if (is.null(maxent4)) maxent4 = maxent1
    maxent <- list(maxent1, maxent2, maxent3, maxent4)
    maxent <- lapply(maxent, raster::values)
    mean_maxent <- sapply(maxent, mean, na.rm=T)
    mean_maxent <- mean_maxent / max(mean_maxent)
    
    if (is.null(abundanceRegion)) {
        abundanceRegion = bathymetry
        abundanceRegion[water_cells] <- 1
    } else {
        abundanceRegion = as.integer(raster::values(abundanceRegion))
    }
    
    if (is.null(fisheryRegion)) {
        fisheryRegion = bathymetry
        fisheryRegion[water_cells] <- 1
    } else {
        fisheryRegion = as.integer(raster::values(fisheryRegion))
    }

    x <- data.table(
        bathymetry = round(bathymetry, precision),
        distToCoast = round(distToCoast, precision),
        abundanceRegion = abundanceRegion,
        fisheryRegion = fisheryRegion,
        block = block,
        foodLevel = foodLevel,
        maxent1 = maxent[[1]],
        maxent2 = maxent[[2]],
        maxent3 = maxent[[3]],
        maxent4 = maxent[[4]])
    
    # replace NAs with integer values, so they can be conveniently handled when reading data in Rcpp.
    x[is.na(bathymetry), bathymetry := 1]
    x[is.na(block), block := -1]
    x[is.na(fisheryRegion), fisheryRegion := -1]
    x[is.na(abundanceRegion), abundanceRegion := -1]
    x[is.na(distToCoast), distToCoast := 0]
    x[is.na(foodLevel), foodLevel := 0]    
    x[is.na(maxent1), maxent1 := 0]
    x[is.na(maxent2), maxent2 := 0]
    x[is.na(maxent3), maxent3 := 0]
    x[is.na(maxent4), maxent4 := 0]

    header = sprintf("%d;%d;%d;%d;%d;%d;%d;%d;%f;%f;%f;%f\n", ncol, nrow, max(x$abundanceRegion), max(x$fisheryRegion), 
                     max(x$block), sum(x$bathymetry<0), as.integer(sqrt(max(table(x$block[x$block!=-1])))), sum(x$foodLevel>0), 
                     mean_maxent[1], mean_maxent[2], mean_maxent[3], mean_maxent[4])
    # write header
    cat(header, file = filename, fill = FALSE, append = FALSE);
    # write data
    data.table::fwrite(x, sep = ";", quote = FALSE, row.names = FALSE, col.names = FALSE,
           file = filename, append = TRUE)
    options(scipen=scipen)
}
