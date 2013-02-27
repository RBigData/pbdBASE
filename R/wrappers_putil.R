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
  
  if (i[2L] < 0){
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


base.ddiagtk <- function(x, descx, reduce=FALSE, proc.dest='all')
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  if (reduce){
    rd <- 'Y'
    
    if (proc.dest[1L] == 'all')
      rdest <- cdest <- -1
    else {
      if (length(proc.dest)==1){
        src <- base.pcoord(ICTXT=descx[2L], PNUM=proc.dest)
        rsrc <- src[[1L]]
        csrc <- src[[2L]]
      }
    }
  }
  else {
    rd <- 'N'
    rdest <- cdest <- -1
  }
  
  ret <- .Call("R_PDDIAGTK", 
               x, as.integer(dim(x)), as.integer(descx), as.integer(min(dx@dim)),
               rd, as.integer(rdest), as.integer(cdest),
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
  
  ret <- new("ddmatrix", Data=out, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
  
  return( ret )
}



base.rcolcpy <- function(dx, dy, xcol, ycol)
{
  ldim <- dx@ldim
  
  descx <- base.descinit(dim=dx@dim, bldim=dx@bldim, ldim=ldim, ICTXT=dx@ICTXT)
  descy <- base.descinit(dim=dy@dim, bldim=dy@bldim, ldim=dy@ldim, ICTXT=dy@ICTXT)
  
  if (!is.double(dx@Data))
    storage.mode(dx@Data) <- "double"
  if (!is.double(dy@Data))
    storage.mode(dy@Data) <- "double"
  
  ret <- .Call("R_RCOLCPY", 
               dx@Data, as.integer(ldim), as.integer(descx), as.integer(xcol), dy@Data, as.integer(descy), as.integer(ycol), as.integer(length(ycol)),
               PACKAGE="pbdBASE")
  
  dx@Data <- ret
  
  return( dx )
}

base.rrowcpy <- function(dx, dy, xrow, yrow)
{
  ldim <- dx@ldim
  
  descx <- base.descinit(dim=dx@dim, bldim=dx@bldim, ldim=ldim, ICTXT=dx@ICTXT)
  descy <- base.descinit(dim=dy@dim, bldim=dy@bldim, ldim=dy@ldim, ICTXT=dy@ICTXT)
  
  if (!is.double(dx@Data))
    storage.mode(dx@Data) <- "double"
  if (!is.double(dy@Data))
    storage.mode(dy@Data) <- "double"
  
  ret <- .Call("R_RROWCPY", 
               dx@Data, as.integer(ldim), as.integer(descx), as.integer(xcol), dy@Data, as.integer(descy), as.integer(ycol), as.integer(length(ycol)),
               PACKAGE="pbdBASE")
  
  dx@Data <- ret
  
  return( dx )
}

