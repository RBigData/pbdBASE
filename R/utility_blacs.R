# Checking if ICTXT is valid
base.valid_context <- function(ICTXT, ..., override=FALSE)
{
  if (!override)
    if (ICTXT==0 || ICTXT==1 || ICTXT==2) 
      comm.stop(paste("Context", ICTXT, "is protected"))
  else if (ICTXT < 0) 
    comm.stop("Negative BLACS context is not allowed")
  else if (as.integer(ICTXT)!=ICTXT) 
    comm.stop("Non-integer BLACS contexts are not permitted")
  else if (!exists(paste(".__blacs_gridinfo_", ICTXT, sep=""))){
    comm.warning(paste("Context", ICTXT, "does not exist"))
    return( invisible(1) )
  } 
}

valid_context <- base.valid_context



# finding the minimum avaliable BLACS context
base.minctxt <- function(after=0)
{
  if (after < -1)
    comm.stop("Error : contexts must be non-negative")
  
  after <- after+1
  
  while(TRUE){
    nm <- paste(".__blacs_gridinfo_", after, sep="")
    if (!exists(nm, envir = .pbdBASEEnv))
      break
    after <- after+1
  }
  
  return(after)
}

minctxt <- base.minctxt



# get BLACS communicator info
base.blacs <- function(ICTXT=0)
{
  ICTXT <- as.integer(ICTXT)
  
  grid <- paste(".__blacs_gridinfo_", ICTXT, sep="")
  
  if (!exists(grid, envir=.pbdBASEEnv))
    comm.stop(paste("Processor grid ICTXT=", ICTXT, " does not exist.  Make sure you called init.grid()", sep=""))
  
  gridinfo <- get(grid, envir=.pbdBASEEnv)
  
  return(gridinfo)
}

blacs <- base.blacs



# ordinary MPI process number from BLACS grid location
base.pnum <- function(ICTXT, PROW, PCOL)
{
  blacs_ <- base.blacs(ICTXT=ICTXT)
#  nprows <- blacs_$NPROW
#  
#  PNUM <- PROW * nprows + PCOL
#  
  NPCOL <- blacs_$NPCOL
  
  PNUM <- PROW * NPCOL + PCOL
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


