base.mksubmat <- function(x, bldim=.BLDIM, ICTXT=0)
{
  if (length(bldim)==1)
    bldim <- rep(bldim, 2)
  
  dim <- dim(x)
  ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT, fixme=TRUE)
  
  descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  subx <- .Call("R_MKSUBMAT", x, as.integer(ldim), as.integer(descx), PACKAGE="pbdBASE")
  
  new("ddmatrix", Data=subx, dim=dim, ldim=ldim, bldim=bldim, CTXT=ICTXT)
}

base.mkgblmat <- function(dx, proc.dest='all')
{
  ICTXT <- dx@CTXT
  
  dim <- dx@dim
  ldim <- dx@ldim
  bldim <- dx@bldim
  
  descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
  
  if (proc.dest=='all')
    rsrc <- csrc <- -1
  else {
    dest <- base.pcoord(ICTXT=ICTXT, PNUM=proc.dest)
    rsrc <- dest[[1]]
    csrc <- dest[[2]]
  }
  
  if (!is.double(dx@Data))
    storage.mode(dx@Data) <- "double"
  
  .Call("R_MKGBLMAT", dx@Data, as.integer(descx), as.integer(rsrc), as.integer(csrc), PACKAGE="pbdBASE")
}


base.tri2zero <- function(dx, uplo='L', diag='N')
{
  ldim <- dx@ldim
  
  descx <- base.descinit(dim=dx@dim, bldim=dx@bldim, ldim=ldim, ICTXT=dx@CTXT)
  
  if (!is.double(dx@Data))
    storage.mode(dx@Data) <- "double"
  
  ret <- .Call("R_PTRI2ZERO", 
               as.character(uplo), as.character(diag), dx@Data, as.integer(ldim), as.integer(descx), 
               PACKAGE="pbdBASE")
  
  dx@Data <- ret
  
  return(dx)
}



