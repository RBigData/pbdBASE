#### Lastest load into a package.

.First.lib <- function(lib, pkg){
  if(! is.loaded("spmd_initialize", PACKAGE = "pbdMPI")){
    library.dynam("pbdMPI", "pbdMPI", lib)
    if(pbdMPI:::comm.is.null(0L) == -1){
      pbdMPI:::init()
    }
  }

  if(! is.loaded("slap_blacs_gridinit", PACKAGE = "pbdSLAP")){
    library.dynam("pbdSLAP", "pbdSLAP", lib)
  }

  library.dynam("pbdBASE", pkg, lib)
} # End of .First.lib().

.Last.lib <- function(libpath){
  ### To free all BLACS points.
  pbdBASE:::finalize(mpi.finalize = FALSE)
  library.dynam.unload("pbdBASE", libpath)
} # End of .Last.lib().

