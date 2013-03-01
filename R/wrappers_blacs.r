# ################################################
# ------------------------------------------------
# Process grid
# ------------------------------------------------
# ################################################

# "Optimal" process grid when nprow and npcol are empty
base.procgrid <- function(nprocs)
{
  out <- .Fortran("OPTIMALGRID", as.integer(nprocs), integer(1), integer(1))
  out[[1L]] <- NULL
  return(list(nprow=out[[2L]], npcol=out[[1L]]))
}

procgrid <- base.procgrid



#Initialize Process Grid
base.init.grid <- function(nprow, npcol, ICTXT)
{
  pbdMPI::init() # initialize pbdMPI communicator
  
  nprocs <- pbdMPI::comm.size()
  
  if (missing(ICTXT)){
    if (exists(".__blacs_gridinfo_0")){
      comm.warning("Context 0 is already initialized. No new grid created")
      return(invisible(1))
    } else {
      ICTXT <- 0
    }
  }
  else if (ICTXT==0 || ICTXT==1 || ICTXT==2) 
    comm.stop("Contexts 0, 1, and 2 are reserved; use 3 or above.")
  else if (ICTXT < 0)
    comm.stop("Context must be at least 3")
  else if (ICTXT - as.integer(ICTXT) > 0) 
    comm.stop("Context must be an integer")
  
  # optimal size grid if parameters are missing
  if (missing(nprow) && missing(npcol)){
    procs <- base.procgrid(nprocs=nprocs)
    nprow <- procs$nprow
    npcol <- procs$npcol
  } 
  else if (missing(nprow) && !missing(npcol))
    comm.stop("You must also provide a value for 'nprow'")
  else if (!missing(nprow) && missing(npcol)) 
    comm.stop("You must also provide a value for 'npcol'")
  else if (nprow*npcol > nprocs) 
    comm.stop(paste("Error: grid size of ", nprow, "*", npcol, " is not possible with ", nprocs, " processes", sep=""))
  
  # Informing the user of creation
  if (ICTXT==0)
    pbdMPI::comm.cat(sprintf("%s", paste("Using ", nprow, "x", npcol, " for the default grid size\n\n", sep="")), quiet=TRUE)
  else
    pbdMPI::comm.cat(sprintf("%s", paste("Grid ICTXT=", ICTXT, " of size ", nprow, "x", npcol, " successfully created\n\n", sep="")), quiet=TRUE)
  
  value <- .Fortran("mpi_blacs_initialize", 
                  NPROW=as.integer(nprow), 
                  NPCOL=as.integer(npcol), 
                  ICTXT=as.integer(0), 
                  MYROW=as.integer(0), 
                  MYCOL=as.integer(0),
                  PACKAGE="pbdBASE"
                )
  
  if (ICTXT==0){
    # Full context
    assign(x=".__blacs_gridinfo_0", value=value, envir=.pbdBASEEnv)
    # Row context
#    if (.__blacs_gridinfo_0$NPROW==1)
#      assign(x=".__blacs_gridinfo_1", 
#       value=value, 
#       envir=.pbdBASEEnv
#      )
#    else
      assign(x=".__blacs_gridinfo_1", 
             value=.Fortran("mpi_blacs_initialize", 
                    NPROW=as.integer(1), NPCOL=as.integer(nprow*npcol), 
                    ICTXT=as.integer(1), MYROW=as.integer(0), 
                    MYCOL=as.integer(0) ),
             envir=.pbdBASEEnv )
    # Col context
      assign(x=".__blacs_gridinfo_2", 
             value=.Fortran("mpi_blacs_initialize", 
                    NPROW=as.integer(nprow*npcol), NPCOL=as.integer(1), 
                    ICTXT=as.integer(2), MYROW=as.integer(0), 
                    MYCOL=as.integer(0) ),
             envir=.pbdBASEEnv )
  }
  else
    assign(x=paste(".__blacs_gridinfo_", ICTXT, sep=""), value=value, envir=.pbdBASEEnv)

  if (!exists(".__blacs_initialized"))
    assign(x=".__blacs_initialized", value=TRUE, envir=.pbdBASEEnv)
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


