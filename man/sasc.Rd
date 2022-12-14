% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sasc.R
\name{sasc}
\alias{sasc}
\title{Creates a serialized (flattened) raster stack from raster files}
\usage{
sasc(
  filename,
  bathymetry,
  distToCoast,
  block = NULL,
  block_size = NULL,
  abundanceRegion = NULL,
  fisheryRegion = NULL,
  foodLevel,
  maxent1,
  maxent2 = NULL,
  maxent3 = NULL,
  maxent4 = NULL,
  precision = 2
)
}
\arguments{
\item{filename}{Output filename}

\item{bathymetry}{Raster. Average depths.}

\item{distToCoast}{Raster. Distance (in units of cells) to the nearest point on the coast}

\item{block}{Raster. Block membership of each cell. Alternatively, use block_size to automatically assign blocks}

\item{block_size}{Single integer. The size of blocks in units of cells. Set this to automatically divide the landscape into blocks. Ignored if \emph{block} is specified}

\item{abundanceRegion}{Raster. Abundance region (if missing, the whole landscape is grouped into one abundance region)}

\item{fisheryRegion}{Raster. fishery region (if missing, the whole landscape is grouped into one fishery region)}

\item{foodLevel}{Raster. Denotes whether a cell is a food patch or not (1 = food patch, 0 otherwise)}
}
\value{
None
}
\description{
Creates a sasc data file from raster files, optionally cropping raster data, that can be used by \code{\link[=sim]{sim()}}
The so-called serialized ASCII (sasc) file is a sort of flattened raster stack, where the cell values
of different rasters are used as column vectors in a data file. The first line in the sasc file is a metadata header
specifying the extent of the landscape, block size, number of food patches and the max values for each of
the different columns (for more convenient processing when reading data again later).
Note that all rasters used must have the same extent and resolution.
}
