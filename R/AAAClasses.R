
#' S4 class to represent IBM simulation results
#' @slot conf List specifying simulation options
#' @slot filename Character specifying filename of serialized ascii data (sasc)
#' @slot data Dataframe containing data from sasc file
#' @slot tracks Dataframe containing movement tracks for any followed porpoises in simulation
#' @slot raster A RasterLayer object based on bathymetric data
setClass("posim", 
         representation(
             conf = "vector",
             filename = "character",
             data = "data.frame",
             bycatch = "data.frame",
             tracks = "data.frame",
             abundance = "data.frame",
             gillnets = "data.frame",
             block = "data.frame",
             structure = "data.frame",
             raster = "RasterLayer"
         )
)

setGeneric("posim", function(steps, ...) {
    standardGeneric("posim", steps = steps, ...)
})

setGeneric("animate", function(x, ...) {
    standardGeneric("animate")
})