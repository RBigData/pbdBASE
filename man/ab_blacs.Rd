\name{BLACS}
\alias{BLACS}
\alias{blacs}
\title{ Get BLACS Context Grid Information }
\description{
  Grabs the existing BLACS context grid information.
}
\usage{
  blacs(ICTXT = 0)
}
\arguments{
  \item{ICTXT}{BLACS context number.}
}
\details{
  BLACS contexts have important internal use, and advanced users familiar with 
  ScaLAPACK might find some advantage in directly manipulating these process 
  grids. Most users should not need to directly manage BLACS contexts, in this 
  function or elsewhere.
  
  The function effectively serves as a shorthand for
  
  \code{eval(parse(text=paste(".__blacs_gridinfo_", ICTXT, sep="")))}
}
\value{
  Returns a list with 5 elements: \code{NPROW} and \code{NPCOL}, the number of
  process rows and columns respectively; \code{ICTXT}, the associated BLACS 
  context number; \code{MYROW} and \code{MYCOL}, the current process' row and 
  column position in the process grid.
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

mygrid <- blacs(0)

comm.print(mygrid)

finalize()
}
}
\keyword{BLACS}
