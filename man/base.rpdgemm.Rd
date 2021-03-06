% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_pblas_level3.r
\name{base.rpdgemm}
\alias{base.rpdgemm}
\title{rpdgemm}
\usage{
base.rpdgemm(transx, transy, x, descx, y, descy, descc)
}
\arguments{
\item{transx, transy}{'T' or 'N' for transpose or not.}

\item{x, y}{Matrix.}

\item{descx, descy, descc}{ScaLAPACK descriptor array.}
}
\description{
Matrix-Matrix Multiply.
}
\details{
For advanced users only. See pbdDMAT for high-level functions.
}
