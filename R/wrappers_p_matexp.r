p_mateye <- function(n, bldim=.BLDIM, ICTXT=.ICTXT)
{
  if (length(bldim==1))
    bldim <- rep(bldim, 2L)
  
  dim <- c(n, n)
  ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT, fixme=TRUE)
  
  desc <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
  
  out <- .Call("R_p_mateye", as.integer(desc), as.integer(ldim), PACKAGE="pbdBASE")
  
  ret <- new("ddmatrix", Data=out, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
}



p_matpow_by_squaring <- function(A, b=1)
{
  b <- as.integer(b)
  
  if (!is.double(A@Data))
    storage.mode(A@Data) <- "double"
  
  desca <- base.descinit(dim=A@dim, bldim=A@bldim, ldim=A@ldim, ICTXT=A@ICTXT)
  
  out <- .Call("R_p_matpow_by_squaring", A@Data, as.integer(desca), b, PACKAGE="pbdBASE")
  
  ret <- new("ddmatrix", Data=out, dim=A@dim, ldim=A@ldim, bldim=A@bldim, ICTXT=A@ICTXT)
  
  return( ret )
}



p_matexp_pade <- function(A)
{
  if (!is.double(A@Data))
    storage.mode(A@Data) <- "double"
  
  desca <- base.descinit(dim=A@dim, bldim=A@bldim, ldim=A@ldim, ICTXT=A@ICTXT)
  
  out <- .Call("R_p_matexp_pade", A@Data, as.integer(desca), PACKAGE="pbdBASE")
  
  N <- new("ddmatrix", Data=out$N, dim=A@dim, ldim=A@ldim, bldim=A@bldim, ICTXT=A@ICTXT)
  D <- new("ddmatrix", Data=out$D, dim=A@dim, ldim=A@ldim, bldim=A@bldim, ICTXT=A@ICTXT)
  
  R <- pbdDMAT::solve(D) %*% N
  
  return( R )
}

