#' matexp
#' 
#' Serial matrix exponentiation.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param A
#' Matrix to exponentiate.
#' @param p
#' Pade' expansion size.
#' @param t
#' Scaling factor.
#' @return exp(A)
#' 
#' @useDynLib pbdBASE R_matexp
#' @export
base.matexp <- function(A, p=6, t=1)
{
  if (!is.double(A))
    storage.mode(A) <- "double"
  
  R <- .Call(R_matexp, A, as.integer(p), as.double(t))
  R
}



#' p_matpow_by_squaring_wrap
#' 
#' Matrix power by squaring.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param A
#' Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' @param b
#' Power.
#' @return A powered matrix.
#' 
#' @useDynLib pbdBASE R_p_matpow_by_squaring
#' @export
base.p_matpow_by_squaring_wrap <- function(A, desca, b=1)
{
  b <- as.integer(b)
  desca <- as.integer(desca)
  
  if (!is.double(A))
    storage.mode(A) <- "double"
  
  ret <- .Call(R_p_matpow_by_squaring, A, desca, b)
  ret
}



#' p_matexp_pade_wrap
#' 
#' Pade' expansion.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param A
#' Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' @param p
#' Order of the Pade' approximation.
#' @return Results of Pade' expansion.
#' 
#' @useDynLib pbdBASE R_p_matexp_pade
#' @export
base.p_matexp_pade_wrap <- function(A, desca, p=6)
{
  desca <- as.integer(desca)
  
  if (!is.double(A))
    storage.mode(A) <- "double"
  
  ret <- .Call(R_p_matexp_pade, A, desca, as.integer(p))
  ret
}
