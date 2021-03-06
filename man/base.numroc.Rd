% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dims.r
\name{base.numroc}
\alias{base.numroc}
\title{numroc}
\usage{
base.numroc(dim, bldim, ICTXT = 0, fixme = TRUE)
}
\arguments{
\item{dim}{Global dim.}

\item{bldim}{Blocking dim.}

\item{ICTXT}{BLACS context.}

\item{fixme}{Should ldims be "rounded" to 0 or not.}
}
\value{
A vector of local dim.
}
\description{
NUMber of Rows Or Columns
}
\details{
For advanced users only. See pbdDMAT for high-level functions.
}
\examples{
spmd.code <- "
  suppressMessages(library(pbdMPI))
  suppressMessages(library(pbdBASE))
  init.grid()

  ### Set desc.
  dim <- c(6L, 5L)
  bldim <- c(3L, 3L)
  ldim <- base.numroc(dim = dim, bldim = bldim)
  comm.print(ldim)

  finalize()
"
pbdMPI::execmpi(spmd.code = spmd.code, nranks = 2L)

}
