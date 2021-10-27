#' rpdlaprnt
#'
#' Matrix printer.
#'
#' For advanced users only. See pbdDMAT for high-level functions.
#'
#' @param m,n
#' Number rows/cols.
#' @param a
#' Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#'
#' @useDynLib pbdBASE R_PDLAPRNT
#' @export
base.rpdlaprnt <- function(m, n, a, desca)
{
  if (!is.double(a))
    storage.mode(a) <- "double"
  
  .Call(R_PDLAPRNT,
        as.integer(m), as.integer(n),
        a, as.integer(desca),
        as.character(deparse(substitute(a))),
        6L #WCC: 0 for stderr, 6 for stdout. Both are disabled.
        )
  
  invisible(0)
}



#' rpdgemr2d
#' 
#' General 2d block cyclic redistribution function.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param x
#' Matrix.
#' @param descx,descy
#' ScaLAPACK descriptor array.
#' 
#' @useDynLib pbdBASE R_PDGEMR2D
#' @export
base.rpdgemr2d <- function(x, descx, descy)
{
  ldimy <- base.numroc(dim=descy[3L:4L], bldim=descy[5L:6L], ICTXT=descy[2L])
  ldimy <- as.integer(ldimy)
  descx <- as.integer(descx)
  descy <- as.integer(descy)
  m <- descx[3L]
  n <- descx[4L]
  
  # context 0 is always passed since pxgemr2d 
  # requires the grids to have at least 1 processor in common
  ### TODO integrate PIGEMR2D
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  ret <- .Call(R_PDGEMR2D, m, n, x, descx, ldimy, descy, 0L)
  
  if (!base.ownany(dim=c(m, n), bldim=descy[5L:6L], ICTXT=descy[2L]))
    ret <- matrix(0.0, 1L, 1L)
  
  ret
}



#' Next Best Divisor
#' 
#' Given integers n and d, with n>d, this function finds the "next
#' best divisor" of n which is greater than or equal to d.
#' 
#' Suprisingly useful for thinking about processor grid shapes.
#' 
#' @param n
#' The divident (number divided into).
#' @param d
#' The candidate divisor.
#' @return The "next best divisor" interger
#' 
#' @examples
#' spmd.code <- "
#'   suppressMessages(library(pbdMPI))
#'   suppressMessages(library(pbdBASE))
#'   init.grid()
#'
#'   base.nbd(100, 10) # 10 divides 100, so 10 is returned
#'   base.nbd(100, 11) # 11 does not, so the 'next best' divisor, 20, is returned
#'
#'   finalize()
#' "
#' pbdMPI::execmpi(spmd.code = spmd.code, nranks = 1L)
#' 
#' @useDynLib pbdBASE R_nbd
#' @export
base.nbd <- function(n, d)
{
  .Call(R_nbd, as.integer(n), as.integer(d))
}



isint <- function(x, epsilon=1e-8)
{
  if (is.numeric(x))
  {
    if (abs(x-as.integer(x)) < epsilon)
      return( TRUE )
    else
      return( FALSE )
  }
  else
    return( FALSE )
}
