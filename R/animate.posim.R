#' Animates simulated movement tracks
#' 
#' Creates an animated GIF file showing in each frame only the most recent movement tracks
#' 
#' @rdname animate
#' @param x The returned object from a call to [posim::posim()]
#' @param select Vector with indices denoting the ids/indices of which individuals to animate. Defaults to NULL, which animates all individuals.
#' @param steps The number of steps to animate.
#' @param tail.length 
#' 
#' @import gifski
#' @md
#' @export
#animate <- function(x, select = 1, steps = NULL, tail.length = 10, blink.duration = 12, delay = 1/24, filename = "track.gif") {
setMethod("animate", signature(x = "posim"), function(x, select = NULL, steps = NULL, tail.length = 10, blink.duration = 12, delay = 1/24, filename = "track.gif") {
    
    # by default, show all steps
    if (is.null(steps)) {
        steps <- 1:max(x@tracks$step)
    } else {
        if (class(steps) != "numeric" | any(steps < 0)) {
            stop("steps is incorrectly specified. See ?posim::animate")
        }
        
        if (length(steps) == 1) steps = 1:steps
    }
    
    z <- x@tracks[x@tracks$step %in% steps & x@tracks$id,]
    z$id <- as.integer(factor(z$id)) # have ids start from 1
    n <- max(z$id)
    
    if (is.null(select)) {
        if (n > 10) {
            warning("follow > 10; only animating 10 random individuals. Use select = 1:", n, " to overwrite");
            select <- sample(z$id, size = 10)
        } else {
            select <- 1:n
        }
    } else if (!all(select %in% unique(z$id))) {
        stop("select must be one or more of 1:", n)
    }
    
    z <- z[z$id %in% select,]
    
    
    #filenames <- sprintf("%s/%s.png", tempdir(), steps)

    filenames <- sprintf("/Volumes/xiaopenyou/Temp/%s.png", steps)
    # graphing setup
    pal <- colorRampPalette(c("lightblue", "orangered1"))
    cols <- pal(200)
    legend_image <- grDevices::as.raster(matrix(cols, nrow=1))
    
    z$col <- cols[1]
    z$col[which(z$energy>0.5)] <- cols[round(z$energy[z$energy>0.5], 1)*10]
    tracks <- split(z, z$id)
    tracks <- lapply(tracks, function(x) {
        found_food <- 1 + which(x$energy[-1] - x$energy[-nrow(x)] > 0)
        max_step <- max(x$step)
        for (step in found_food) {
            x$col[seq(step, min(step+blink.duration, max_step))] <- "purple"
        }
        x
    })
    
    nstep <- length(steps)
    ntracks <- length(tracks)
    cex.tail <- seq(0.45, 0.225, -0.025)
    tail <- integer()
    
    for (i in steps) {
        message(sprintf("\rProcessing step %d/%d...", i, max(steps)), appendLF = FALSE)
        if (i > 1) {
            tail <- seq(i-1, max(1, i-tail.length), -1)
        }
    
        png(filename = filenames[i], width = 800, height = 800)
        raster::plot(x@raster, main = sprintf("sim day %.02d (step %d)", 1+floor(steps[i]/48), steps[i]))
        usr <- par()$usr
        
        for (j in 1:ntracks) {
            if (i <= nrow(tracks[[j]])) {
                points(tracks[[j]]$x[i], tracks[[j]]$y[i], pch=19, cex = 0.5, col = tracks[[j]]$col[i])
                if (length(tail) > 0) {
                    points(tracks[[j]]$x[tail], tracks[[j]]$y[tail], pch=19, cex = cex.tail[1:length(tail)], col = tracks[[j]]$col[tail])
                }
            }
        }
        
        rasterImage(legend_image, usr[1] + (usr[2]-usr[1])*0.05, 
                    usr[3] + (usr[4]-usr[3])*0.93, 
                    usr[1] + (usr[2]-usr[1])*0.20, 
                    usr[3] + (usr[4]-usr[3])*0.95, 0)
        
        text(usr[1] + (usr[2]-usr[1])*0.05, usr[3] + (usr[4]-usr[3])*0.96, "Energy level", cex = 0.5, adj =0)
        text(usr[1] + (usr[2]-usr[1])*0.04, usr[3] + (usr[4]-usr[3])*0.94, "0", cex = 0.5, adj=0.5)
        text(usr[1] + (usr[2]-usr[1])*0.21, usr[3] + (usr[4]-usr[3])*0.94, "20", cex = 0.5, adj=0.5)
        
        dev.off()
    }
    message("\rCreating GIF...", i, max(steps), appendLF = FALSE)
    
    gifski::gifski(filenames, gif_file = filename, delay = 1/24, width=800, height = 1000)
    unlink(filenames)
    message(sprintf("\rAll done! GIF saved as %s", filename), appendLF = FALSE)
})