base.mksubmat <- function(x, bldim=.BLDIM, ICTXT=0)
{
  if (length(bldim)==1)
    bldim <- rep(bldim, 2)
  
  dim <- dim(x)
  ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT, fixme=TRUE)
  
  descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  subx <- .Call("R_MKSUBMAT", 
                x, as.integer(ldim), as.integer(descx), 
                PACKAGE="pbdBASE")
  
  new("ddmatrix", Data=subx, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
}


base.mkgblmat <- function(dx, proc.dest='all')
{
  ICTXT <- dx@ICTXT
  
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
  
  ret <- .Call("R_MKGBLMAT", 
       dx@Data, as.integer(descx), as.integer(rsrc), as.integer(csrc), 
       PACKAGE="pbdBASE")
  
  return( ret )
  
}



base.dallreduce <- function(dx, op='sum', scope='All')
{
  ldim <- dx@ldim
  
  descx <- base.descinit(dim=dx@dim, bldim=dx@bldim, ldim=ldim, ICTXT=dx@ICTXT)
  
  if (!is.double(dx@Data))
    storage.mode(dx@Data) <- "double"
  
  .Call("R_DALLREDUCE", 
        dx@Data, as.integer(ldim), as.integer(descx), as.character(op), as.character(scope),
        PACKAGE = 'pbdBASE')
}


base.tri2zero <- function(dx, uplo='L', diag='N')
{
  ldim <- dx@ldim
  
  descx <- base.descinit(dim=dx@dim, bldim=dx@bldim, ldim=ldim, ICTXT=dx@ICTXT)
  
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
  
  descx <- base.descinit(dim=dx@dim, bldim=dx@bldim, ldim=ldim, ICTXT=dx@ICTXT)
  
  if (!is.double(dx@Data))
    storage.mode(dx@Data) <- "double"
  
  if (!is.double(vec))
    storage.mode(vec) <- "double"
  
  ret <- .Call("R_PDSWEEP", 
               dx@Data, as.integer(ldim), as.integer(descx), vec, as.integer(length(vec)), as.integer(MARGIN), as.character(FUN),
               PACKAGE="pbdBASE")
  
  dx@Data <- ret
  
  return(dx)
}


base.rl2blas <- function(dx, vec, FUN)
{
  ldim <- dx@ldim
  
  descx <- base.descinit(dim=dx@dim, bldim=dx@bldim, ldim=ldim, ICTXT=dx@ICTXT)
  
  if (!is.double(dx@Data))
    storage.mode(dx@Data) <- "double"
  
  if (!is.double(vec))
    storage.mode(vec) <- "double"
  
  ret <- .Call("R_RL2BLAS", 
               dx@Data, as.integer(ldim), as.integer(descx), vec, as.integer(length(vec)), as.integer(FUN),
               PACKAGE="pbdBASE")
  
#####  dx@Data <- ret
#####  
#####  return(dx)
  return(ret)
}


base.rl2insert <- function(dx, vec, i, j)
{
  ldim <- dx@ldim
  
  descx <- base.descinit(dim=dx@dim, bldim=dx@bldim, ldim=ldim, ICTXT=dx@ICTXT)
  
  if (all(i<0)){
    new <- 1:dx@dim[1]
    i <- new[-which(new %in% abs(i))] # FIXME make this less stupid
  }
  
  if (all(j<0)){
    new <- 1:dx@dim[2]
    j <- new[-which(new %in% abs(j))] # FIXME make this less stupid
  }
  
  if (!is.double(dx@Data))
    storage.mode(dx@Data) <- "double"
  
  if (!is.double(vec))
    storage.mode(vec) <- "double"
  
  ret <- .Call("R_RL2INSERT", 
               dx@Data, as.integer(ldim), as.integer(descx), vec, as.integer(length(vec)), as.integer(i), as.integer(length(i)), as.integer(j), as.integer(length(j)),
               PACKAGE="pbdBASE")
  
  dx@Data <- ret
  
  return(dx)
}


base.ddiagtk <- function(dx)
{
  ldim <- dx@ldim
  
  descx <- base.descinit(dim=dx@dim, bldim=dx@bldim, ldim=ldim, ICTXT=dx@ICTXT)
  
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

