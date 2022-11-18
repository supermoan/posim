#' Plot simulation results
#' 
#' Plots simulation results, either simulated abundance per month, or simulated movement tracks for one or more individuals. 
#' The background map is created with a call to [raster::plot()] on the bathymetry data. For the time being, landscape coordinates
#' are local to the grid size used, and are not translated back into their real-world representation.
#' 
#' @rdname plot
#' @export
#' @param x The returned object from a call to [posim::posim()].
#' @param y Ignored. 
#' @param tracks Logical. If true, plots simulated movement tracks.
#' @param select Vector with indices denoting the ids of individuals to plot. Use NULL (the default) to plot up to 10 random tracks.
#' @param highlight If specified, ids that match will be plotted in red, while all other tracks are colored black.
#' @param ... Further graphical parameters from par, passed to [raster::plot()].
#' 
#' @details
#' Note that the plotting function allocates new ids to all porpoises by calling as.integer(factor(x@tracks$id)).
#' 
#' @return None
#'
#' @seealso [posim::sim()], [posim::animate()], [posim::sasc()]
#' @examples 
#' \dontrun{
#' # plot change in total abundance over time
#' plot(x)
#'
#' # plot tracks for all porpoises
#' plot(x, tracks = TRUE)
#'
#' # plot the first 3 tracks
#' plot(x, select = 1:3)
#'
#' # plot the first 5 tracks, and highlight numbers 1 and 2.
#' plot(x, select = 1:5, highlight = c(1, 2))
#' }
#' 
#' @md

setMethod("plot", signature(x = "posim", y = "character"), function(x, y) {
    
    p <- par(no.readonly = TRUE)
    par(mar = c(5.1,4.2,4.1,2.1))
    a <- x@abundance
    s <- x@structure
    
    rf <- max(a$N)/max(a$food, na.rm=T)
    re <- max(a$N)/20
    
    plot(a$N, type="l", ylim=c(0, max(a$N)*1.25), ylab="\nValue", xlab="Month", las = 1,
         main = "Simulated abundance per month")
    at <- pretty(c(0, max(a$N)) / re) 
    axis(side = 4, at = at*re, labels = at, las=1)
    mtext(side = 3, line = 0.5, cex = 0.75,
          sprintf("foodGrowthRate = %.02f, distEnergyMultiplier = %.02f, maxU = %.02f", 
                  x@conf$foodGrowthRate, x@conf$distEnergyMultiplier, x@conf$maxU))
    abline(h = pretty(1:max(a$N)), v = pretty(1:length(a$N)), col = "grey80")
    
    # add polygons for each ageClass
    xpos <- 1:nrow(a)
    lasty <- rep(0, nrow(a))
    cols <- c("grey60", "grey40", "black")
    
    for (i in seq(0, max(s$ageClass))) {
        y <- lasty + s$N[s$ageClass == i]
        lines(xpos, y, col = cols[i+1], lwd = 2)
        lasty <- y
    }
    
    lines(a$N, col = "black", lwd=4)
    lines(a$food * rf, col = "red", lwd=2)
    lines(a$energy * re, col = "blue", lwd=2)
    legend("topright", inset = 0.05, lty = 1, lwd=2, 
           col = c(rev(cols), "red", "blue"), 
           bg = "transparent",
           legend = c("N (all ages)", "N (adults + juveniles; age < 10)", "N (juveniles; age < 4)", sprintf("%.03f x food", rf), sprintf("%.02f x average energy", re)))
    
    par(p)
    return(invisible(NULL))
})

setMethod("plot", signature(x = "posim", y = "missing"), function(x, tracks = TRUE, steps = NULL, select = NULL, highlight = NULL, legend = FALSE, ...) {
    
    if ((missing(tracks) & missing(select) & missing(highlight)) | tracks == FALSE) {
        plot(x, y = "abundance")
        return(invisible(NULL))
    }
    
    if (nrow(x@tracks) == 0) {
        stop("No tracks to plot! Rerun simulation with follow enabled. See ?posim")
    }
    
    if (is.null(steps)) {
        steps <- 1:max(x@tracks$step)
    } else {
        if (class(steps) != "integer" | any(steps < 0)) {
            stop("steps is incorrectly specified. See ?posim::plot")
        }
        
        if (length(steps) == 1) steps = 1:steps
    }
    
    
    x@tracks <- x@tracks[x@tracks$step %in% steps,]
    x@tracks$id <- as.integer(factor(x@tracks$id))
    n <- max(x@tracks$id)
    
    if (is.null(select)) {
        if (n > 10) {
            warning("follow > 10; only plotting 10 random individuals. Use select = 1:", n, " to overwrite");
            select <- sample(x@tracks$id, size = 10)
        } else {
            select <- 1:n
        }
    } else if (!all(select %in% unique(x@tracks$id))) {
        stop("select must be one or more of 1:", n)
    }

    tracks <- x@tracks[x@tracks$id %in% select,]

    main <- sprintf("Simulated tracks for %d individuals (step %d-%d)", length(select), min(steps), max(steps))
    raster::plot(x@raster, maxpixels = raster::ncell(x@raster), las = 1, main = main, ...)
    
    if (x@conf$steps == 0) { # no steps run - only plotting starting locations
        points(tracks$x, tracks$y, pch = 4, ...)
    } else {
        # split tracks by porpoise id
        tracks <- split(tracks, tracks$id)
        
        jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                         "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
        
        if (!is.null(highlight)) {
            if (all(highlight %in% names(tracks))) {
                cols <- rep("black", times = length(tracks))
                cols[highlight] <- "red"
            } else {
                stop("highlight must be one or more of ", min(select), ":", max(select))
            }
        } else {
            cols <- jet.colors(length(tracks))
        }
        
        # make sure highlighted tracks are on top
        j <- 1:length(tracks)
        j <- j[!j %in% highlight]
        j <- c(j, highlight)
        
        for (i in j) {
            n <- nrow(tracks[[i]]) # how many moves did this porp make?
            if (n > 2) {
                lines(tracks[[i]]$x, tracks[[i]]$y, col = cols[i])
            }
            #points(track$x[n], track$y[n], pch=19, cex = 0.1, col = cols[unique(track$id)])
        }
    }
    
    legend("topleft", legend = sprintf("porp%d", sort(j)), col = cols, lty = 1, inset = 0.05, cex=0.8)
    return(invisible(NULL))
})
