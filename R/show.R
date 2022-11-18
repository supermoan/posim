setMethod("show", signature(object = "posim"), function(object) {
    
    dN = tail(object@abundance$N, 1) / sum(object@conf$N)
    year = substr(as.character(Sys.Date()),1,4)
    start <- as.Date(object@conf$start -1, origin = sprintf("%s-01-01", year))
    stop <- as.Date(object@conf$start -1 + floor(object@conf$steps/48), origin = sprintf("%s-01-01", year))
    
    growth <- aggregate(cbind(N = object@abundance$N[-1]), 
                                    by = list(year = 1 + ceiling(object@abundance$step[-1]/(48*365))),
                                    FUN = tail, n = 1)
    
    growth <- rbind(cbind(year = 0, N = object@abundance$N[1]), growth)
    
    printf("class:     : %s\n", class(object))
    with(object@conf$header, {
        printf("extent     : %d, %d, %d (nrow, ncol, ncell)\n", ymax, xmax, xmax * ymax)
        printf("landscape  : %d, %d (traversable cells, food patches)\n", ntrav, nfood)
        printf("regions    : %d, %d (abundance regions, fishery regions)\n", nabund, nfish)
        printf("blocks     : %d blocks (each %d x %d)\n", nblock, block_size, block_size)
        printf("maxent     : %.03f, %.03f, %.03f, %.03f\n", maxent1, maxent2, maxent3, maxent4)
    })
    
    printf("sim time   : %s, %s (from, to)\n", start, stop)
    printf("sim length : %d, %.02f (30-minute steps, days)\n", object@conf$steps, round(object@conf$steps/48,2))
    printf("tracks     : %d\n", length(object@conf$follow))
    printf("abundance  : %d, %d, %s%.01f%% (start, stop, %% change)\n", sum(object@conf$N), tail(object@abundance$N, 1), ifelse(dN > 1, "+", ""), (dN-1)*100)
    printf("rates      : %.02f%% (average yearly growth rate)\n", mean((growth$N[-1]-growth$N[-nrow(growth)]) / growth$N[-nrow(growth)]))
    printf("filename   : %s (%.01f Mb)\n", object@filename, file.size(object@filename)/1024/1024)
    printf("fishdata   : %s\n", ifelse(!is.null(object@conf$fish), object@conf$fish, NA))
    
    #printf(" First day of simulation: day %d (%s)", object@conf$start, as.Date(object@conf$start, origin = sprintf("%d-01-01", substr(as.character(Sys.Date()),1,4)))
          
})