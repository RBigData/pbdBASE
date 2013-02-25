# R version of ScaLAPACK tool DESCINIT
# Creates ScaLAPACK descriptor array for distributed matrix
base.descinit <- function(dim, bldim, ldim, ICTXT=0)
{
  desc <- integer(9)
  desc[1] <- 1L       # matrix type --- 1 for dense
  desc[2] <- ICTXT    # CTXT_A
  desc[3] <- dim[1]   # M_A
  desc[4] <- dim[2]   # N_A
  desc[5] <- bldim[1] # MB_A
  desc[6] <- bldim[2] # NB_A
  desc[7] <- 0L       # RSRC_A
  desc[8] <- 0L       # CSRC_A
  desc[9] <- ldim[1]  # LLD_A
  
  return(desc)
}

# R version of ScaLAPACK tool NUMROC
# Computes ldim of ddmatrix given dim and bldim
# if fixme=TRUE, then returned local dimensions which are < 1 will 
# be corrected (made = 1), and otherwise will be left unchanged.
base.numroc <- function(dim, bldim, ICTXT=0, fixme=TRUE)
{
  
  blacs_ <- base.blacs(ICTXT=ICTXT)
  
  MYP <- c(blacs_$MYROW, blacs_$MYCOL)
  PROCS <- c(blacs_$NPROW, blacs_$NPCOL)

  ISRCPROC <- 0
  
  ldim <- numeric(2)
  for (i in 1:2){
    MYDIST <- (PROCS[i] + MYP[i] -  ISRCPROC) %% PROCS[i]
    NBLOCKS <- floor(dim[i] / bldim[i])
    ldim[i] <- floor(NBLOCKS / PROCS[i]) * bldim[i]
    EXTRABLKS <- NBLOCKS %% PROCS[i]

    if (is.na(EXTRABLKS))
      EXTRABLKS <- 0

    if (MYDIST < EXTRABLKS)
      ldim[i] <- ldim[i] + bldim[i]
    else if (MYDIST == EXTRABLKS)
      ldim[i] <- ldim[i] + dim[i] %% bldim[i]
  }

  if (fixme){
    if (any(is.na(ldim)))
      ldim[which(is.na(ldim))] <- 0
    if (any(ldim<1)) ldim <- c(1,1) # FIXME
  }

  return(ldim)
}

numroc <- base.numroc


# For use with local arithmetic; basically does nothing if the 
# local storage is just filler to make scalapack happy
# return is logical answer to the question:  'do I own anything?', 
# and not a C-style return
base.ownany <- function(dim, bldim, ICTXT=0)
{
  if (length(bldim)==1)
    bldim <- rep(bldim, 2)
  check <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT, fixme=FALSE)
  
  if (any(check<1))
    return(FALSE)
  else
    return(TRUE)
}

ownany <- function(x, ..., dim, bldim, ICTXT=0)
{
  if (!missing(bldim) && length(bldim)==1)
    bldim <- rep(bldim, 2)
  if (!missing(x) && is.ddmatrix(x))
    return( base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT) )
  else if (!missing(dim) && !missing(bldim) && missing(x) && is.numeric(dim) && is.numeric(bldim))
    return( base.ownany(dim=dim, bldim=bldim, ICTXT=ICTXT) )
  else{
    print("Error: bad input(s) in ownany()")
    stop("")
  }
}


# Hook into ScaLAPACK tool PDLAPRNT
base.rpdlaprnt <- function(m, n, a, desca)
{
  if (!is.double(a))
    storage.mode(a) <- "double"
  
  .Call("R_PDLAPRNT", 
        as.integer(m), as.integer(n),
        dx@Data, as.integer(desca),
        as.character(deparse(substitute(dx))),
        6L,  #WCC: 0 for stderr, 6 for stdout. Both are disabled.
        PACKAGE="pbdBASE"
        )
  
  return( invisible(0) )
}

