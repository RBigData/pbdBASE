#' ScaLAPACK Wrappers and Utilities
#' 
#' A package contains the basic methods for dealing with distributed data
#' types, as well as the data types themselves.
#' 
#' \tabular{ll}{ Package: \tab pbdBASE\cr Type: \tab Package\cr License: \tab
#' MPL\cr LazyLoad: \tab yes\cr } This package requires an MPI library
#' (OpenMPI, MPICH2, or LAM/MPI).
#' 
#' @import utils methods pbdSLAP
#' @importFrom pbdMPI allreduce comm.print comm.stop comm.rank comm.warning comm.is.null bcast
#' 
#' @name pbdBASE-package
#' @docType package
#' @author Drew Schmidt \email{wrathematics AT gmail.com}, Wei-Chen Chen, George
#' Ostrouchov, and Pragneshkumar Patel.
#' @references Programming with Big Data in R Website: \url{https://pbdr.org/}
#' @keywords Package
NULL
