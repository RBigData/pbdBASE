# "Optimal" process grid when nprow and npcol are empty
base.procgrid <- function(nprocs)
{
  out <- .Call("optimal_process_grid", as.integer(nprocs), PACKAGE="pbdBASE")
#  return(list(nprow=out[1], npcol=out[2]))
  return(list(nprow=out[2], npcol=out[1]))
}

procgrid <- base.procgrid


base.valid_context <- function(ICTXT, ..., override=FALSE)
{
  if (!override)
    if (ICTXT==0 || ICTXT==1 || ICTXT==2) {
      comm.print(paste("Context", ICTXT, "is protected"))
      stop("")
     }
  else if (ICTXT < 0) {
    comm.print("Negative BLACS context is not allowed")
    stop("")
  }
  else if (as.integer(ICTXT)!=ICTXT) {
    comm.print("Non-integer BLACS contexts are not permitted")
    stop("")
  }
  else if (!exists(paste(".__blacs_gridinfo_", ICTXT, sep=""))){
    warning(paste("Context", ICTXT, "does not exist"))
    return( invisible(1) )
  } 
}

valid_context <- base.valid_context


base.minctxt <- function(after=0)
{
  if (after < 0){
    comm.print("Error : contexts must be positive")
    stop("")
  }
  
  after <- after+1
  
  while(TRUE){
    x <- try( blacs(after), silent=T )
    if (inherits(x=x, what="try-error"))
      break
    after <- after+1
  }
  
  return(after)
}

minctxt <- base.minctxt

#Initialize Process Grid
base.init.grid <- function(nprow, npcol, ICTXT)
{
  pbdMPI::init() # initialize pbdMPI communicator
  
  nprocs <- pbdMPI::comm.size()

  if (missing(ICTXT)){
    if (exists(".__blacs_gridinfo_0")){
      warning("Context 0 is already initialized. No new grid created")
      return(invisible(1))
    } else {
      ICTXT <- 0
    }
  } else if (ICTXT==0 || ICTXT==1 || ICTXT==2) {
    comm.print("Contexts 0, 1, and 2 are reserved; use 3 or above.")
    stop("")
  } else if (ICTXT < 0){
    comm.print("Context must be at least 3")
    stop("")
  } else if (ICTXT - as.integer(ICTXT) > 0) {
    comm.print("Context must be an integer")
    stop("")
  }

  # optimal size grid if parameters are missing
  if (missing(nprow) && missing(npcol)){
    procs <- base.procgrid(nprocs=nprocs)
    nprow <- procs$nprow
    npcol <- procs$npcol
  } else if (missing(nprow) && !missing(npcol)) {
    comm.print("You must also provide a value for 'nprow'")
    stop("")
  }
  else if (!missing(nprow) && missing(npcol)) {
    comm.print("You must also provide a value for 'npcol'")
    stop("")
  }
  else if (nprow*npcol > nprocs) {
    comm.print(paste("Error: grid size of ", nprow, "*", npcol, " is not possible with ", nprocs, " processes", sep=""))
    stop("")
  }

  # Informing the user of creation
  if (ICTXT==0)
    pbdMPI::comm.cat(sprintf("%s", paste("Using a default grid size of ", nprow, "x", npcol, "\n\n", sep="")), quiet=TRUE)
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
    assign(x=".__blacs_gridinfo_0", value=value, envir=.GlobalEnv)
    # Row context
#    if (.__blacs_gridinfo_0$NPROW==1)
#      assign(x=".__blacs_gridinfo_1", 
#       value=value, 
#       envir=.GlobalEnv
#      )
#    else
      assign(x=".__blacs_gridinfo_1", 
             value=.Fortran("mpi_blacs_initialize", 
                    NPROW=as.integer(1), NPCOL=as.integer(nprow*npcol), 
                    ICTXT=as.integer(1), MYROW=as.integer(0), 
                    MYCOL=as.integer(0) ),
             envir=.GlobalEnv )
    # Col context
#    if (.__blacs_gridinfo_0$NPCOL==1)
#      assign(x=".__blacs_gridinfo_2", 
#       value=value, 
#       envir=.GlobalEnv
#      )
#    else
      assign(x=".__blacs_gridinfo_2", 
             value=.Fortran("mpi_blacs_initialize", 
                    NPROW=as.integer(nprow*npcol), NPCOL=as.integer(1), 
                    ICTXT=as.integer(2), MYROW=as.integer(0), 
                    MYCOL=as.integer(0) ),
             envir=.GlobalEnv )
  }
  else
    assign(x=paste(".__blacs_gridinfo_", ICTXT, sep=""), value=value, envir=.GlobalEnv)

  if (!exists(".__blacs_initialized"))
    assign(x=".__blacs_initialized", value=TRUE, envir=.GlobalEnv)
  invisible(0) # quiet return
}

init.grid <- base.init.grid


# get BLACS communicator info
base.blacs <- function(ICTXT=0)
{
  ICTXT <- as.integer(ICTXT)
  gridinfo <- get(paste(".__blacs_gridinfo_", ICTXT, sep=""), envir=.GlobalEnv)
  
  return(gridinfo)
}

blacs <- base.blacs


# shut down a BLACS context
base.gridexit <- function(ICTXT, ..., override=FALSE)
{
  base.valid_context(ICTXT=ICTXT, override=override)
  
  blacs_ <- base.blacs(ICTXT=ICTXT)
  FCTXT <- blacs_$ICTXT

  if (blacs_$MYROW != -1 && blacs_$MYCOL != -1)
    .Fortran("BLACS_GRIDEXIT", ICONTXT=as.integer(FCTXT), PACKAGE="pbdBASE")

  rm(list = paste(".__blacs_gridinfo_", ICTXT, sep=""), envir=.GlobalEnv)

  return( invisible(0) )
}

gridexit <- base.gridexit


# ordinary MPI process number from BLACS grid location
base.pnum <- function(ICTXT, PROW, PCOL)
{
  blacs_ <- base.blacs(ICTXT=ICTXT)
  nprows <- blacs_$NPROW
  
  PNUM <- PROW * nprows + PCOL
  
  return( PNUM )
}

pnum <- base.pnum


# BLACS grid location from ordinary MPI process number
base.pcoord <- function(ICTXT, PNUM)
{
  blacs_ <- blacs(ICTXT=ICTXT)
  nprows <- blacs_$NPROW
  
  PROW <- as.integer(PNUM / nprows)
  PCOL <- PNUM %% nprows
  
  return( list(PROW=PROW, PCOL=PCOL) )
}

pcoord <- base.pcoord


# row/column sums of matrices via BLACS
base.blacs.sum <- function(SCOPE, A, dim, na.rm=FALSE, ICTXT=0, means=FALSE, num=1) # SCOPE= 'Row', 'Column', 'All'
{
  if (SCOPE==1)
    SCOPE <- 'Row'
  if (SCOPE==2)
    SCOPE <- 'Column'

  if (SCOPE=='Row')
    if (!means)
      f <- function(x, na.rm=na.rm) rowSums(x, na.rm=na.rm)
    else
      f <- function(x, na.rm=na.rm) rowSums(x, na.rm=na.rm) / num
  else if (SCOPE=='Column')
    if (!means)
      f <- function(x, na.rm=na.rm) colSums(x, na.rm=na.rm)
    else
      f <- function(x, na.rm=na.rm) colSums(x, na.rm=na.rm) / num

  A <- f(x=A, na.rm=na.rm)

  M <- max(dim)
  N <- 1
  LDA <- length(A)
  
  mxm <- pbdMPI::allreduce(M, op='max')
  if (length(A)==1 && A[1]==0){
    A <- numeric(mxm)
    M <- mxm
  }

  ### WCC: out should be allocated within .Call.
  out <- .Call("R_row_col_sums",
               as.integer(ICTXT), as.character(SCOPE),
               as.integer(M), as.integer(N), A, as.integer(LDA),
               PACKAGE="pbdBASE"
               )

  return(out)
}

blacs.sum <- base.blacs.sum

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
  if (exists(".__blacs_initialized", envir = .GlobalEnv)){
    base.blacsexit(CONT=TRUE)
    rm(list = ".__blacs_initialized", envir = .GlobalEnv)
  }
    
  pbdMPI::finalize(mpi.finalize=mpi.finalize)
}

