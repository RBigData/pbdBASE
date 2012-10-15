\name{Gridexit}
\alias{Gridexit}
\alias{gridexit}
\title{ gridexit }
\description{
  Frees a BLACS context.
}
\usage{
  gridexit(ICTXT, ..., override = FALSE)
}
\arguments{
  \item{ICTXT}{BLACS context number.}
  \item{...}{additional arguments.}
  \item{override}{logical; if TRUE, ignores normal check preventing
  the closing of \code{ICTXT} values of 0, 1, and 2.}
}
\details{
  For advanced users only. 
  
  The function frees the requested BLACS context. It is a trivial
  wrapper for the BLACS routine \code{BLACS_GRIDEXIT}. Also removes
  the globally stored R object \code{.__blacs_gridinfo_ICTXT}. 
  
  Contexts 0, 1, and 2 can not be freed in this way unless the
  argument \code{override=FALSE}. This will probably break something
  and I do not recommend it.
}
\value{
  Silently returns 0 when successful. Silently returns 1 when 
  requested \code{ICTXT} does not exist.
}
\seealso{
  \code{\link{InitGrid}}
}
\keyword{BLACS}
\keyword{Global Variables}
