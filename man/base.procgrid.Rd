% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base_mpi_blacs.r
\name{base.procgrid}
\alias{base.procgrid}
\title{procgrid}
\usage{
base.procgrid(nprocs)
}
\arguments{
\item{nprocs}{Number of processors.}
}
\value{
A list contains nprow and npcol.
}
\description{
"Optimal" process grid when nprow and npcol are empty
}
\details{
For advanced users only. See pbdDMAT for high-level functions.
}
\examples{
spmd.code <- "
  suppressMessages(library(pbdMPI))
  suppressMessages(library(pbdBASE))
  init.grid()

  opt <- base.procgrid(4)
  comm.print(opt)

  opt <- base.procgrid(6)
  comm.print(opt)

  opt <- base.procgrid(8)
  comm.print(opt)

  finalize()
"
pbdMPI::execmpi(spmd.code = spmd.code, nranks = 1L)

}
