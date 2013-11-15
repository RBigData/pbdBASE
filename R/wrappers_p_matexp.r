p_mateye <- function(n, bldim=.BLDIM, ICTXT=.ICTXT)
{
  if (length(bldim==1))
    bldim <- rep(bldim, 2L)
  
  dim <- c(n, n)
  ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT, fixme=TRUE)
  
  desc <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
  
  out <- .Call("R_p_mateye", desc, PACKAGE="pbdDMAT")
  
  ret <- new("ddmatrix", Data=out, dim=dim, ldim=dim(Data), bldim=bldim, ICTXT=ICTXT)
}
