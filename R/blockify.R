
blockify <- function(r, block_size) {
    nx <- seq(1, nrow(r), block_size) # first column number in each block
    ny <- seq(1, ncol(r), block_size) # first row number in each block
    block <- 1 # block number, incremented in each block
    
    for (x in nx) {
        for (y in ny) {
            xsec <- seq(x, min(x+block_size-1, nrow(r)))
            ysec <- seq(y, min(y+block_size-1, ncol(r)))
            r[xsec, ysec] <- block
            block <- block +1
        }
    }
    r
}