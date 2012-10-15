\name{InitGrid}
\alias{InitGrid}
\alias{init.grid}
\title{ Initialize Process Grid }
\description{
  Manages the creation of BLACS context grids.
}
\usage{
  init.grid(nprow, npcol, ICTXT)
}
\arguments{
  \item{nprow}{number of process rows.}
  \item{npcol}{number of process columns.}
  \item{ICTXT}{BLACS context number.}
}
\details{
  If \code{ICTXT} is missing, three global variables will be created:
  
  .__blacs_gridinfo_0
  
  .__blacs_gridinfo_1
  
  .__blacs_gridinfo_2

  These global variables store the BLACS process grid information for the BLACS
  context corresponding to the trailing digit of the variable. Most users should
  invoke \code{init.grid()} in this fashion, namely with ICTXT missing, and only
  do so once.
  
  Contexts 0, 1, 2, and 3 are reserved. Additional custom contexts
  are possible to create, but they must be integers > 3.
  
  Context 0 is the ``full'' process grid of \code{nprow} by \code{npcol} 
  processes; contexts 1 is the process grid consisting of 1 process row and 
  \code{nprow}*\code{npcol} processes columns; context 2 is the process grid
  consisting of \code{nprow}*\code{npcol} processes rows and 1 process column. 
  These contexts can be redundant depending on the number of prcesses available.
  
  BLACS contexts have important internal use, and advanced users familiar with 
  ScaLAPACK might find some advantage in directly manipulating these process 
  grids. Most users should not need to directly manage BLACS contexts, in this 
  function or elsewhere.
  
  If the \code{nprow} and \code{npcol} values are missing, then a best process
  grid will be chosen for the user based on the total available number of
  processes. Here ``best'' means as close to a square grid as possible.
  
  The global variables \code{.__blacs_gridinfo_ICTXT} are just storage mechanisms
  to avoid needing to directly invoke the BLACS routine \code{BLACS_GRIDINFO}.

  Additionally, another global variable is created when running 
  \code{init.grid()}, namely \code{.__blacs_initialized}. Its existence 
  is to alert \code{finalize()} to shut down BLACS communicators, if
  necessary, to prevent memory leaks.
}
\value{
  Silently returns 0 when successful. Additionally, several global variables
  are produced.  See Details section.
}
\seealso{
  \code{\link{as.ddmatrix}, \link{BLACS}, \link{Distribute}}
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run with 2 processors by
# > mpiexec -np 2 Rscript demo.r

library(pbdBASE, quiet = TRUE)
init.grid()

finalize()
}
}
\keyword{BLACS}
\keyword{Global Variables}
