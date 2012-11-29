\name{BLACS Exit}
\alias{BLACS Exit}
\alias{blacsexit}
\title{ BLACS Exit }
\description{
  Shuts down all BLACS communicators.
}
\usage{
  blacsexit(CONT = TRUE)
}
\arguments{
  \item{CONT}{logical; determines whether or not to shut down 
  \emph{all} MPI communicators}
}
\details{
  If the user wishes to shut down BLACS communicators but still have 
  access to MPI, then call this function with \code{CONT=TRUE}. 
  Calling \code{blacsexit(CONT=FALSE)} will shut down all MPI 
  communicators, equivalent to calling 
  
  \code{> blacsexit(CONT=TRUE)}
  \code{> finalize(mpi.finalize=TRUE)}
  
  This function is automatically invoked if BLACS communicators are 
  running and \code{finalize()} is called.
}
\value{
  Has an invisible return of 0 when successful.
}
\seealso{
  \code{\link{InitGrid}}
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run with 2 processors by
# > mpiexec -np 2 Rscript demo.r

library(pbdBASE, quiet = TRUE)
init.grid()

blacsexit()

finalize()
}
}
\keyword{BLACS}
