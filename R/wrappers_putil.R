base.mksubmat <- function(x, descx)
{
  ldim <- base.numroc(dim=descx[3L:4L], bldim=descx[5L:6L], ICTXT=descx[2L], fixme=TRUE)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  subx <- .Call("R_MKSUBMAT", 
                x, as.integer(ldim), as.integer(descx), 
                PACKAGE="pbdBASE")
  
  return( subx )
}


base.mkgblmat <- function(x, descx, rsrc, csrc)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  ret <- .Call("R_MKGBLMAT", 
       x, as.integer(descx), as.integer(rsrc), as.integer(csrc), 
       PACKAGE="pbdBASE")
  
  return( ret )
  
}


base.dallreduce <- function(x, descx, op='sum', scope='All')
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  ret <- .Call("R_DALLREDUCE", 
        x, as.integer(dim(x)), as.integer(descx), as.character(op), as.character(scope),
        PACKAGE = 'pbdBASE')
  
  return( ret )
}


base.tri2zero <- function(x, descx, uplo='L', diag='N')
{
  uplo <- toupper(uplo)
  diag <- toupper(diag)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  ret <- .Call("R_PTRI2ZERO", 
               uplo, diag, x, as.integer(dim(x)), as.integer(descx), 
               PACKAGE="pbdBASE")
  
  return( ret )
}


base.pdsweep <- function(x, descx, vec, MARGIN, FUN)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  if (!is.double(vec))
    storage.mode(vec) <- "double"
  
  ret <- .Call("R_PDSWEEP", 
               x, as.integer(dim(x)), as.integer(descx), vec, as.integer(length(vec)), as.integer(MARGIN), as.character(FUN),
               PACKAGE="pbdBASE")
  
  return( ret )
}


base.rl2blas <- function(x, descx, vec, FUN)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  if (!is.double(vec))
    storage.mode(vec) <- "double"
  
  ret <- .Call("R_RL2BLAS", 
               x, as.integer(dim(x)), as.integer(descx), vec, as.integer(length(vec)), as.integer(FUN),
               PACKAGE="pbdBASE")
  
  return(ret)
}

# matrix-vector insertion
base.rl2insert <- function(x, descx, vec, i, j)
{
  dim <- descx[3L:4L]
  
  if (i[1L] < 0){
    new <- 1L:dim[1L]
    i <- new[-which(new %in% abs(i))] # FIXME make this less stupid
  }
  
  if (j[1L] < 0){
    new <- 1L:dim[2L]
    j <- new[-which(new %in% abs(j))] # FIXME make this less stupid
  }
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  if (!is.double(vec))
    storage.mode(vec) <- "double"
  
  ret <- .Call("R_RL2INSERT", 
               x, as.integer(dim(x)), as.integer(descx), vec, as.integer(length(vec)), as.integer(i), as.integer(length(i)), as.integer(j), as.integer(length(j)),
               PACKAGE="pbdBASE")
  
  return( ret )
}


base.ddiagtk <- function(x, descx, proc.dest='all')
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  if (proc.dest[1L] == 'all')
    rdest <- cdest <- -1
  else {
    if (length(proc.dest)==1){
      src <- base.pcoord(ICTXT=descx[2L], PNUM=proc.dest)
      rsrc <- src[[1L]]
      csrc <- src[[2L]]
    }
  }
  
  ldiag <- min(descx[3L:4L])
  
  ret <- .Call("R_PDGDGTK", 
               x, as.integer(dim(x)), as.integer(descx), as.integer(ldiag),
               as.integer(rdest), as.integer(cdest),
               PACKAGE="pbdBASE")
  
  return( ret )
}


base.ddiagmk <- function(diag, descx)
{
  ldim <- base.numroc(dim=descx[3L:4L], bldim=descx[5L:6L], ICTXT=descx[2L])
  
  if (!is.double(diag))
    storage.mode(diag) <- "double"
  
  out <- .Call("R_PDDIAGMK", 
               as.integer(ldim), as.integer(descx), diag, as.integer(length(diag)),
               PACKAGE="pbdBASE")
  
  return( out )
}



base.rcolcpy <- function(x, descx, y, descy, xcol, ycol)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  ret <- .Call("R_RCOLCPY", 
               x, as.integer(dim(x)), as.integer(descx), as.integer(xcol), y, as.integer(descy), as.integer(ycol), as.integer(length(ycol)),
               PACKAGE="pbdBASE")
  
  return( ret )
}

base.rcolcpy2 <- function(x, descx, y, descy, xcol, ycol)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  ret <- .Call("R_RCOLCPY2", 
               x, as.integer(dim(x)), as.integer(descx), as.integer(xcol), as.integer(length(xcol)), y, as.integer(descy), as.integer(ycol), as.integer(length(ycol)),
               PACKAGE="pbdBASE")
  
  return( ret )
}

base.rrowcpy <- function(x, descx, y, descy, xrow, yrow)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  ret <- .Call("R_RROWCPY", 
               x, as.integer(dim(x)), as.integer(descx), as.integer(xrow), y, as.integer(descy), as.integer(yrow), as.integer(length(yrow)),
               PACKAGE="pbdBASE")
  
  return( ret )
}

base.rrowcpy2 <- function(x, descx, y, descy, xrow, yrow)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  ret <- .Call("R_RROWCPY2", 
               x, as.integer(dim(x)), as.integer(descx), as.integer(xrow), as.integer(length(xrow)), y, as.integer(descy), as.integer(yrow), as.integer(length(yrow)),
               PACKAGE="pbdBASE")
  
  return( ret )
}

