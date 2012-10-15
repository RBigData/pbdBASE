\name{Coords}
\alias{Coords}
\alias{pnum}
\alias{pcoord}
\title{ Interchange Between Process Number and BLACS Coordinates }
\description{
  Grabs the existing BLACS context grid information.
}
\usage{
  pnum(ICTXT, PROW, PCOL)
  pcoord(ICTXT, PNUM)
}
\arguments{
  \item{ICTXT}{BLACS context number.}
  \item{PROW, PCOL}{BLACS grid location row/column}
  \item{PNUM}{process rank}
}
\details{
  For advanced users only. These functions are simple recreations of
  the BLACS routines \code{BLACS_PNUM} and \code{BLACS_PCOORD}. The
  former gets the process number associated with the BLACS process 
  grid location \code{c(MYPROW, MYPCOL)}, while the latter does the 
  reverse.
}
\value{
  \code{pnum} returns an integer; \code{pcoord} returns a list 
  containing elements \code{PROW} and \code{PCOL}.
}
\seealso{
  \code{\link{BLACS}, \link{InitGrid}}
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run with 2 processors by
# > mpiexec -np 2 Rscript demo.r

library(pbdBASE, quiet = TRUE)
init.grid()

blacs_ <- blacs(ICTXT = 0)

# get the ICTXT = 0 BLACS coordsinates for process 0
myCoords <- pcoord(ICTXT = 0, PNUM = 0)

comm.print(myCoords)

finalize()
}
}
\keyword{BLACS}
\keyword{Global Variables}
