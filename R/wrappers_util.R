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
  
  if (proc.dest[1]=='all')
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



base.dallreduce <- function(dx, op='sum', scope='All')
{
  ldim <- dx@ldim
  
  descx <- base.descinit(dim=dx@dim, bldim=dx@bldim, ldim=ldim, ICTXT=dx@CTXT)
  
  if (!is.double(dx@Data))
    storage.mode(dx@Data) <- "double"
  
  .Call("R_DALLREDUCE", 
        dx@Data, as.integer(ldim), as.integer(descx), as.character(op), as.character(scope),
        PACKAGE = 'pbdBASE')
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


base.pdsweep <- function(dx, vec, MARGIN, FUN)
{
  ldim <- dx@ldim
  
  descx <- base.descinit(dim=dx@dim, bldim=dx@bldim, ldim=ldim, ICTXT=dx@CTXT)
  
  if (!is.double(dx@Data))
    storage.mode(dx@Data) <- "double"
  
  ret <- .Call("R_PTRI2ZERO", 
               dx@Data, as.integer(ldim), as.integer(descx), as.double(vec), as.integer(length(vec)), as.integer(MARGIN), as.character(FUN),
               PACKAGE="pbdBASE")
  
  dx@Data <- ret
  
  return(dx)
}


base.ddiagtk <- function(dx)
{
  ldim <- dx@ldim
  
  descx <- base.descinit(dim=dx@dim, bldim=dx@bldim, ldim=ldim, ICTXT=dx@CTXT)
  
  if (!is.double(dx@Data))
    storage.mode(dx@Data) <- "double"
  
  ret <- .Call("R_PDDIAGTK", 
               dx@Data, as.integer(ldim), as.integer(descx), as.integer(min(dx@dim)),
               PACKAGE="pbdBASE")
  
  
  return( ret )
}


base.ddiagmk <- function(diag, nrow, ncol, bldim, ICTXT=0)
{
  dim <- c(nrow, ncol)
  ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
  
  descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
  
  if (!is.double(diag))
    storage.mode(diag) <- "double"
  
  out <- .Call("R_PDDIAGMK", 
               as.integer(ldim), as.integer(descx), diag, as.integer(length(diag)),
               PACKAGE="pbdBASE")
  
  ret <- new("ddmatrix", Data=out, dim=dim, ldim=ldim, bldim=bldim, CTXT=ICTXT)
  
  return( ret )
}



