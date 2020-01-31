#' (Un)Distribute
#' 
#' (Un)Distribute matrix.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param x
#' Matrix.
#' @param descx
#' ScaLAPACK descriptor array.
#' 
#' @examples
#' spmd.code <- "
#'   suppressMessages(library(pbdMPI))
#'   suppressMessages(library(pbdBASE))
#'   init.grid()
#'
#'   ### Set data matrix and desc.
#'   x <- matrix(as.double(1:30), nrow = 6, ncol = 5)
#'   dim <- dim(x)
#'   bldim <- c(3L, 3L)
#'   ldim <- base.numroc(dim = dim, bldim = bldim)
#'   descx <- base.descinit(dim = dim, bldim = bldim, ldim = ldim)
#'
#'   ### Redistribute from rank 0.
#'   dx <- base.mksubmat(x, descx)
#'   comm.print(dx, all.rank = TRUE)
#'
#'   finalize()
#' "
#' pbdMPI::execmpi(spmd.code = spmd.code, nranks = 2L)
#' 
#' @useDynLib pbdBASE R_MKSUBMAT
#' @rdname lclgblmat
#' @export
base.mksubmat <- function(x, descx)
{
  ldim <- base.numroc(dim=descx[3L:4L], bldim=descx[5L:6L], ICTXT=descx[2L], fixme=TRUE)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  subx <- .Call(R_MKSUBMAT, x, as.integer(ldim), as.integer(descx))
  subx
}



#' @param rsrc,csrc
#' Row/column source.
#' 
#' @useDynLib pbdBASE R_MKGBLMAT
#' @rdname lclgblmat
#' @export
base.mkgblmat <- function(x, descx, rsrc, csrc)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  ret <- .Call(R_MKGBLMAT, 
       x, as.integer(descx), as.integer(rsrc), as.integer(csrc))
  
  ret
}



#' tri2zero
#' 
#' Zero Triangle
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param x
#' Matrix.
#' @param descx
#' ScaLAPACK descriptor array.
#' @param uplo
#' Triangle.
#' @param diag
#' Zero diagonal as well.
#' 
#' @useDynLib pbdBASE R_PTRI2ZERO
#' @export
base.tri2zero <- function(x, descx, uplo='L', diag='N')
{
  uplo <- toupper(uplo)
  diag <- toupper(diag)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  ret <- .Call(R_PTRI2ZERO, 
               uplo, diag, x, as.integer(dim(x)), as.integer(descx))
  
  ret
}



#' pdsweep
#' 
#' Matrix-Vector Sweep
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param x
#' Matrix.
#' @param descx
#' ScaLAPACK descriptor array.
#' @param vec
#' Vector
#' @param MARGIN
#' Rows or columns.
#' @param FUN
#' Function.
#' 
#' @useDynLib pbdBASE R_PDSWEEP
#' @export
base.pdsweep <- function(x, descx, vec, MARGIN, FUN)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  if (!is.double(vec))
    storage.mode(vec) <- "double"
  
  ret <- .Call(R_PDSWEEP, 
               x, as.integer(dim(x)), as.integer(descx), vec, as.integer(length(vec)), as.integer(MARGIN), as.character(FUN))
  
  ret
}



#' diag
#' 
#' Grab diagonal or create distributed diagonal matrix.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param x
#' Matrix.
#' @param descx
#' ScaLAPACK descriptor array.
#' @param proc.dest
#' Who owns the result.
#' @return diagonal elements of matrix or a diagonal matrix
#' 
#' @examples
#' spmd.code <- "
#'   suppressMessages(library(pbdMPI))
#'   suppressMessages(library(pbdBASE))
#'   init.grid()
#'
#'   ### Set data matrix and desc.
#'   x <- matrix(as.double(1:25), nrow = 5, ncol = 5)
#'   dim <- dim(x)
#'   bldim <- c(3L, 3L)
#'   ldim <- base.numroc(dim = dim, bldim = bldim)
#'   descx <- base.descinit(dim = dim, bldim = bldim, ldim = ldim)
#'
#'   ### Get diagonal
#'   diag.x <- base.ddiagtk(x, descx)
#'   comm.print(diag.x)
#'
#'   finalize()
#' "
#' pbdMPI::execmpi(spmd.code = spmd.code, nranks = 2L)
#'
#' @useDynLib pbdBASE R_PDGDGTK
#' @name diag
#' @rdname diag
#' @export
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
  
  ret <- .Call(R_PDGDGTK, 
               x, as.integer(dim(x)), as.integer(descx), as.integer(ldiag),
               as.integer(rdest), as.integer(cdest))
  
  ret
}

#' @param diag
#' Diagonal.
#' 
#' @examples
#' spmd.code <- "
#'   suppressMessages(library(pbdMPI))
#'   suppressMessages(library(pbdBASE))
#'   init.grid()
#'
#'   ### Set data matrix and desc.
#'   x <- matrix(as.double(1:25), nrow = 5, ncol = 5)
#'   dim <- dim(x)
#'   bldim <- c(3L, 3L)
#'   ldim <- base.numroc(dim = dim, bldim = bldim)
#'   descx <- base.descinit(dim = dim, bldim = bldim, ldim = ldim)
#'
#'   ### Set diagonal
#'   diag.x <- base.ddiagtk(x, descx)
#'   new.x <- base.ddiagmk(diag.x, descx)
#'   comm.print(new.x, all.rank = TRUE)
#'
#'   finalize()
#' "
#' pbdMPI::execmpi(spmd.code = spmd.code, nranks = 2L)
#'
#' @useDynLib pbdBASE R_PDDIAGMK
#' @rdname diag
#' @export
base.ddiagmk <- function(diag, descx)
{
  ldim <- base.numroc(dim=descx[3L:4L], bldim=descx[5L:6L], ICTXT=descx[2L])
  
  if (!is.double(diag))
    storage.mode(diag) <- "double"
  
  out <- .Call(R_PDDIAGMK, 
               as.integer(ldim), as.integer(descx), diag, as.integer(length(diag)))
  
  out
}



#' dhilbmk
#' 
#' Create Hilbert matrix.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param n
#' Size.
#' 
#' @useDynLib pbdBASE R_DHILBMK
#' @export
base.dhilbmk <- function(n)
{
  n <- as.integer(n)
  
  ret <- .Call(R_DHILBMK, n)
  ret
}



#' pdhilbmk
#' 
#' Create Hilbert matrix.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param descx
#' ScaLAPACK descriptor matrix.
#' 
#' @useDynLib pbdBASE R_PDHILBMK
#' @export
base.pdhilbmk <- function(descx)
{
  descx <- as.integer(descx)
  ldim <- as.integer(base.numroc(dim=descx[3L:4L], bldim=descx[5L:6L], ICTXT=descx[2L], fixme=TRUE))
  
  ret <- .Call(R_PDHILBMK, ldim, descx)
  ret
}



#' pdmkcpn1
#' 
#' Create Companion Matrix
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param coef
#' Coefficients vector.
#' @param descx
#' ScaLAPACK descriptor array.
#' 
#' @useDynLib pbdBASE R_PDMKCPN1
#' @export
base.pdmkcpn1 <- function(coef, descx)
{
  ldim <- base.numroc(dim=descx[3L:4L], bldim=descx[5L:6L], ICTXT=descx[2L])
  
  if (!is.double(coef))
    storage.mode(coef) <- "double"
  
  out <- .Call(R_PDMKCPN1, as.integer(ldim), as.integer(descx), coef)
  out
}
