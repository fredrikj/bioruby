neumann <- function(obs,densityfile) {
    densitymatrix <- as.matrix(read.table(densityfile))
    w <- diag(obs) %*% densitymatrix
    w <- w / sum(diag(w))
    e <- eigen(w, only.values=TRUE)$values
    - sum(e * log(e),na.rm=TRUE)
}
#m <- array(c(2,-2,-2,6),dim=c(2,2))
#o <- c(0.8,0.2)
#neumann(o)

