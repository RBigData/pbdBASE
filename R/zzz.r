#### Lastest load into a package.

### Export Namespace does not use .First.lib() and .Last.lib(), but use
### .onLoad() and .onUnload().
# .First.lib <- function(lib, pkg){
# } # End of .First.lib().

.Last.lib <- function(libpath){
  ### To free all BLACS points.
  pbdBASE:::finalize(mpi.finalize = FALSE)
} # End of .Last.lib().

.onLoad <- function(libname, pkgname){
  if(! is.loaded("spmd_initialize", PACKAGE = "pbdMPI")){
    library.dynam("pbdMPI", "pbdMPI", libname)
    if(pbdMPI:::comm.is.null(0L) == -1){
      pbdMPI:::init()
    }
  }

  if(! is.loaded("slap_blacs_gridinit", PACKAGE = "pbdSLAP")){
    library.dynam("pbdSLAP", "pbdSLAP", libname)
  }

  library.dynam("pbdBASE", pkgname, libname)
  invisible()
} # End of .onLoad().

.onUnload <- function(libpath){
  library.dynam.unload("pbdBASE", libpath)
  invisible()
} # End of .onUnload().

