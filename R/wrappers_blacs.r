# ################################################
# ------------------------------------------------
# Process grid
# ------------------------------------------------
# ################################################

# "Optimal" process grid when nprow and npcol are empty
base.procgrid <- function(nprocs)
{
#  out <- .Fortran("OPTIMALGRID", as.integer(nprocs), integer(1), integer(1))
#  out[[1L]] <- NULL
#  return(list(nprow=out[[2L]], npcol=out[[1L]]))
  .Call("R_optimal_grid", as.integer(nprocs), PACKAGE="pbdBASE")
}

procgrid <- base.procgrid



isint <- function(x){
  if (is.numeric(x)){
    if (x-as.integer(x) == 0)
      return( TRUE )
    else
      return( FALSE )
  }
  else
    return( FALSE )
}


# Initialize Process Grid --- these functions have side effects and no return
base.blacs_gridinit <- function(ICTXT, NPROW, NPCOL, ..., quiet = FALSE)
{
  if (missing(ICTXT))
    ICTXT <- base.minctxt(after=-1)
  else if (!isint(x=ICTXT) || ICTXT < 0)
    comm.stop("ICTXT must be a non-negative integer")
  
  
#  pbdMPI::init() # initialize pbdMPI communicator
  nprocs <- pbdMPI::comm.size()
  
  if (missing(NPROW) && missing(NPCOL)){
    procs <- base.procgrid(nprocs=nprocs)
    NPROW <- as.integer(procs$nprow)
    NPCOL <- as.integer(procs$npcol)
  } 
  else if (missing(NPROW) && !missing(NPCOL))
    comm.stop("You must also provide a value for 'NPROW'")
  else if (!missing(NPROW) && missing(NPCOL)) 
    comm.stop("You must also provide a value for 'NPCOL'")
  else if (!isint(x=NPROW) || !isint(x=NPCOL))
    comm.stop("'NPROW' and 'NPCOL' must be integers")
  else if (NPROW*NPCOL > nprocs) 
    comm.stop(paste("Error: grid size of ", NPROW, "x", NPCOL, " is not possible with ", nprocs, " processes", sep=""))
  
  
  nm <- paste(".__blacs_gridinfo_", ICTXT, sep="")
  
  if (exists(nm, envir=.pbdBASEEnv)){
    comm.warning(paste("Context", ICTXT, "is already in use. No new grid created"))
    return( invisible(1) )
  }
  
  value <- .Call("R_blacs_init", 
                 as.integer(NPROW), as.integer(NPCOL), as.integer(ICTXT),
                 PACKAGE = "pbdBASE")
  
  assign(x=nm, value=value, envir=.pbdBASEEnv )
  
  if (ICTXT==0 && !quiet)
    pbdMPI::comm.cat(sprintf("%s", paste("Using ", NPROW, "x", NPCOL, " for the default grid size\n\n", sep="")), quiet=TRUE)
  else if (ICTXT > 0 && !quiet)
    pbdMPI::comm.cat(sprintf("%s", paste("Grid ICTXT=", ICTXT, " of size ", NPROW, "x", NPCOL, " successfully created\n", sep="")), quiet=TRUE)
  
  if (!exists(".__blacs_initialized", envir=.pbdBASEEnv))
    assign(x=".__blacs_initialized", value=TRUE, envir=.pbdBASEEnv)
  
  invisible( 0 )
}

blacs_gridinit <- base.blacs_gridinit



base.init.grid <- function(NPROW, NPCOL, ICTXT, ..., quiet = FALSE)
{
  # initialize pbdMPI communicator
  pbdMPI::init() 
  
  # determine the ICTXT if it is missing
  if (missing(ICTXT)){
    ICTXT <- base.minctxt(after=-1)
  } else if (ICTXT==0 || ICTXT==1 || ICTXT==2){
      comm.stop("Contexts 0, 1, and 2 are reserved; use 3 or above.")
  }
  
  # determine number processor rows/columns
  if (missing(NPROW) && missing(NPCOL)){
    if (exists(".__blacs_gridinfo_0")){
      comm.warning("Context 0 is already initialized. No new grid created")
      return( invisible(1) )
    }
  } else if (missing(NPROW) || missing(NPCOL)){
    comm.stop("You must supply either both 'NPROW' and 'NPCOL' or neither.")
  }
  
  # initialize grid
  base.blacs_gridinit(ICTXT=ICTXT, NPROW=NPROW, NPCOL=NPCOL, quiet=quiet)
  
  
  if (ICTXT==0){
    if (missing(NPROW) && missing(NPCOL)){
      nprocs <- pbdMPI::comm.size()
      
      procs <- base.procgrid(nprocs=nprocs)
      NPROW <- as.integer(procs$nprow)
      NPCOL <- as.integer(procs$npcol)
    } 
    base.blacs_gridinit(ICTXT=1L, NPROW=1, NPCOL=NPROW*NPCOL, quiet=TRUE)
    base.blacs_gridinit(ICTXT=2L, NPROW=NPROW*NPCOL, NPCOL=1, quiet=TRUE)
  }
  
  invisible(0) # quiet return
}

init.grid <- base.init.grid



# shut down a BLACS context
base.gridexit <- function(ICTXT, ..., override=FALSE)
{
  base.valid_context(ICTXT=ICTXT, override=override)
  
  blacs_ <- base.blacs(ICTXT=ICTXT)
  FCTXT <- blacs_$ICTXT

  if (blacs_$MYROW != -1 && blacs_$MYCOL != -1)
    .Fortran("BLACS_GRIDEXIT", ICONTXT=as.integer(FCTXT), PACKAGE="pbdBASE")

  rm(list = paste(".__blacs_gridinfo_", ICTXT, sep=""), envir=.pbdBASEEnv)

  return( invisible(0) )
}

gridexit <- base.gridexit



# exit the blacs grid
base.blacsexit <- function(CONT=TRUE)
{
  .Fortran("BLACS_EXIT", as.integer(CONT), PACKAGE="pbdBASE")
  
  return( invisible(0) )
}

blacsexit <- base.blacsexit



# replacement for pbdMPI::finalize() that automatically shuts BLACS down
finalize <- function(mpi.finalize=.SPMD.CT$mpi.finalize)
{
  if (exists(".__blacs_initialized", envir = .pbdBASEEnv)){
    base.blacsexit(CONT=TRUE)
    rm(list = ".__blacs_initialized", envir = .pbdBASEEnv)
  }
    
  pbdMPI::finalize(mpi.finalize=mpi.finalize)
}


# ################################################
# ------------------------------------------------
# Reductions et al
# ------------------------------------------------
# ################################################

# Sums
base.igsum2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.integer(x))
    storage.mode(x) <- "integer"
  
  out <- .Call("R_igsum2d1", as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST), 
                PACKAGE="pbdBASE")
  
  return( out )
}

base.dgsum2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call("R_dgsum2d1", as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST), 
                PACKAGE="pbdBASE")
  
  return( out )
}


# Max value
base.igamx2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.integer(x))
    storage.mode(x) <- "integer"
  
  out <- .Call("R_igamx2d1", as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST), 
                PACKAGE="pbdBASE")
  
  return( out )
}

base.dgamx2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call("R_dgamx2d1", as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST), 
                PACKAGE="pbdBASE")
  
  return( out )
}


# Min value
base.igamn2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.integer(x))
    storage.mode(x) <- "integer"
  
  out <- .Call("R_igamn2d1", as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST), 
                PACKAGE="pbdBASE")
  
  return( out )
}

base.dgamn2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call("R_dgamn2d1", as.integer(ICTXT), as.character(SCOPE), 
                as.integer(m), as.integer(n), x, as.integer(lda), 
                as.integer(RDEST), as.integer(CDEST), 
                PACKAGE="pbdBASE")
  
  return( out )
}


# point to point communication
base.dgesd2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call("R_dgesd2d1", as.integer(ICTXT), as.integer(m), as.integer(n), 
                x, as.integer(lda), as.integer(RDEST), as.integer(CDEST), 
                PACKAGE="pbdBASE")
  
  return( out )
}

base.dgerv2d <- function(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
{
  if (!is.matrix(x) && !is.vector(x))
    comm.stop("ERROR : object 'x' must be of class matrix or vector.")
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  out <- .Call("R_dgerv2d1", as.integer(ICTXT), as.integer(m), as.integer(n), 
                x, as.integer(lda), as.integer(RDEST), as.integer(CDEST), 
                PACKAGE="pbdBASE")
  
  return( out )
}


