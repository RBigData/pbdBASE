### S4 methods for R copycats
setGeneric(name="sweep", useAsDefault=sweep)
setGeneric(name="print", useAsDefault=print)
setGeneric(name="nrow", useAsDefault=nrow)
setGeneric(name="ncol", useAsDefault=ncol)
setGeneric(name="as.matrix", useAsDefault=as.matrix)
setGeneric(name="na.exclude", useAsDefault=na.exclude)
setGeneric(name="all.equal", useAsDefault=all.equal)

setGeneric(name="as.vector", 
  function(x, ...)
    standardGeneric("as.vector"),
  package="pbdBASE"
)

setGeneric(name="rbind", 
  function(..., ICTXT=0, deparse.level=1)
    standardGeneric("rbind"),
  package="pbdBASE"
)

setGeneric(name="cbind", 
  function(..., ICTXT=0, deparse.level=1)
    standardGeneric("cbind"),
  package="pbdBASE"
)

### S4 methods for new things
setGeneric(name="as.ddmatrix", 
  function(x, ...) 
    standardGeneric("as.ddmatrix"), 
  package="pbdBASE"
)

setGeneric(name="submatrix", 
  function(x, ...) 
    standardGeneric("submatrix"), 
  package="pbdBASE"
)

setGeneric("submatrix<-", 
  function(x, value)
    standardGeneric("submatrix<-"),
  package="pbdBASE"
)

setGeneric(name="ldim", 
  function(x, ...) 
    standardGeneric("ldim"), 
  package="pbdBASE"
)

setGeneric(name="bldim", 
  function(x, ...) 
    standardGeneric("bldim"), 
  package="pbdBASE"
)

setGeneric(name="ctxt", 
  function(x, ...) 
    standardGeneric("ctxt"), 
  package="pbdBASE"
)
