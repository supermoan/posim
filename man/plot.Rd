% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot,posim,character-method}
\alias{plot,posim,character-method}
\title{Plot simulation results}
\usage{
\S4method{plot}{posim,character}(x, y)
}
\arguments{
\item{x}{The returned object from a call to \code{\link[=posim]{posim()}}.}

\item{y}{Ignored.}

\item{tracks}{Logical. If true, plots simulated movement tracks.}

\item{select}{Vector with indices denoting the ids of individuals to plot. Use NULL (the default) to plot up to 10 random tracks.}

\item{highlight}{If specified, ids that match will be plotted in red, while all other tracks are colored black.}

\item{...}{Further graphical parameters from par, passed to \code{\link[raster:plot]{raster::plot()}}.}
}
\value{
None
}
\description{
Plots simulation results, either simulated abundance per month, or simulated movement tracks for one or more individuals.
The background map is created with a call to \code{\link[raster:plot]{raster::plot()}} on the bathymetry data. For the time being, landscape coordinates
are local to the grid size used, and are not translated back into their real-world representation.
}
\details{
Note that the plotting function allocates new ids to all porpoises by calling as.integer(factor(x@tracks$id)).
}
\examples{
\dontrun{
# plot change in total abundance over time
plot(x)

# plot tracks for all porpoises
plot(x, tracks = TRUE)

# plot the first 3 tracks
plot(x, select = 1:3)

# plot the first 5 tracks, and highlight numbers 1 and 2.
plot(x, select = 1:5, highlight = c(1, 2))
}

}
\seealso{
\code{\link[=sim]{sim()}}, \code{\link[=animate]{animate()}}, \code{\link[=sasc]{sasc()}}
}
