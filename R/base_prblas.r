#' Level 2 R-like BLAS
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param x
#' Matrix.
#' @param descx
#' ScaLAPACK descriptor array.
#' @param vec
#' Global vector.
#' @param FUN
#' Function.
#' 
#' @useDynLib pbdBASE R_RL2BLAS
#' @export
base.rl2blas <- function(x, descx, vec, FUN)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  if (!is.double(vec))
    storage.mode(vec) <- "double"
  
  ret <- .Call(R_RL2BLAS, 
               x, as.integer(dim(x)), as.integer(descx), vec, as.integer(length(vec)), as.integer(FUN))
  
  ret
}



#' R-like Matrix-Vector Insertion
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param x
#' Matrix.
#' @param descx
#' ScaLAPACK descriptor array.
#' @param vec
#' Global vector.
#' @param i,j
#' Indices.
#' 
#' @useDynLib pbdBASE R_RL2INSERT
#' @export
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
  
  ret <- .Call(R_RL2INSERT, 
               x, as.integer(dim(x)), as.integer(descx), vec, as.integer(length(vec)), as.integer(i), as.integer(length(i)), as.integer(j), as.integer(length(j)))
  
  ret
}



#' R Column Copy
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param x,y
#' Matrix.
#' @param descx,descy
#' ScaLAPACK descriptor array.
#' @param xcol,ycol
#' Columns.
#' 
#' @useDynLib pbdBASE R_RCOLCPY
#' @export
base.rcolcpy <- function(x, descx, y, descy, xcol, ycol)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  ret <- .Call(R_RCOLCPY, 
               x, as.integer(dim(x)), as.integer(descx), as.integer(xcol), y, as.integer(descy), as.integer(ycol), as.integer(length(ycol)))
  
  ret
}



#' R Column Copy-2
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param x,y
#' Matrix.
#' @param descx,descy
#' ScaLAPACK descriptor array.
#' @param xcol,ycol
#' Columns.
#' 
#' @useDynLib pbdBASE R_RCOLCPY2
#' @export
base.rcolcpy2 <- function(x, descx, y, descy, xcol, ycol)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  ret <- .Call(R_RCOLCPY2, 
               x, as.integer(dim(x)), as.integer(descx), as.integer(xcol), as.integer(length(xcol)), y, as.integer(descy), as.integer(ycol), as.integer(length(ycol)))
  
  ret
}



#' R Row Copy
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param x,y
#' Matrix.
#' @param descx,descy
#' ScaLAPACK descriptor array.
#' @param xrow,yrow
#' Rows.
#' 
#' @useDynLib pbdBASE R_RROWCPY
#' @export
base.rrowcpy <- function(x, descx, y, descy, xrow, yrow)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  ret <- .Call(R_RROWCPY, 
               x, as.integer(dim(x)), as.integer(descx), as.integer(xrow), y, as.integer(descy), as.integer(yrow), as.integer(length(yrow)))
  
  ret
}



#' R Row Copy-2
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param x,y
#' Matrix.
#' @param descx,descy
#' ScaLAPACK descriptor array.
#' @param xrow,yrow
#' Rows.
#' 
#' @useDynLib pbdBASE R_RROWCPY2
#' @export
base.rrowcpy2 <- function(x, descx, y, descy, xrow, yrow)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  ret <- .Call(R_RROWCPY2, 
               x, as.integer(dim(x)), as.integer(descx), as.integer(xrow), as.integer(length(xrow)), y, as.integer(descy), as.integer(yrow), as.integer(length(yrow)))
  
  ret
}



#' R-like Matrix-Vector Sum
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param x
#' Matrix.
#' @param y
#' Vector.
#' @param descx,descy
#' ScaLAPACK descriptor array.
#' 
#' @useDynLib pbdBASE R_PDMVSUM
#' @export
base.pdmvsum <- function(x, descx, y, descy)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  ret <- .Call(R_PDMVSUM, 
               x, as.integer(dim(x)), as.integer(descx), y, as.integer(descy))
  
  ret
}
