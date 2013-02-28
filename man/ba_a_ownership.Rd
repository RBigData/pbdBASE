\name{OwnerShip}
\alias{Ownership}

\alias{base.ownany}

\alias{numroc}

\title{ Determining Local Ownership of a Distributed Matrix }
\description{
  aa
}
\usage{
  base.ownany(dim, bldim, ICTXT = 0)
  numroc(dim, bldim, ICTXT = 0, fixme = TRUE)
}
\arguments{
  \item{dim}{global dimension}
  \item{bldim}{blocking dimension}
  \item{ICTXT}{BLACS context}
  \item{fixme}{logical, controls correction of local dimension return}
}
\details{
  For advanced users only.
  
  \code{numroc()} is a re-implementation at the R level of the 
  ScaLAPACK subroutine NUMROC, which returns the local dimension of 
  the matrix storage, i.e. the dimension for the \code{Data} slot 
  of the distributed matrix on that process. The \code{fixme=} option,
  if \code{TRUE}, returns a minimum of 1 for each dimension. If 
  \code{fixme=FALSE}, then values less than 1 for either dimension 
  are possible, and in this case indicate a lack of local ownership 
  of the global matrix.
  
  \code{ownany()} is a simple wrapper of numroc. The return is the 
  answer to the question 'do I own any of the global matrix?'. 
  Passing a distributed matrix is allowed, but often it is convenient
  to determine that information without even having a distributed 
  matrix on hand. In this case, explicitly passing the appropriate
  information to the arguments \code{dim=}, \code{bldim=} (and 
  \code{CTXT=} as necessary, since it defaults to 0) while leaving
  \code{x} missing will produce the desired result. See the examples
  below for more clarity.
  
  The return for each function is local.
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

iown <- ownany(dim=c(4, 4), bldim=c(2, 2), CTXT=0)
comm.print(iown, all.rank=T)

finalize()
}
}
\keyword{BLACS}
\keyword{Distributing Data}
