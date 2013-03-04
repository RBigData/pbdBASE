# R version of ScaLAPACK tool DESCINIT
# Creates ScaLAPACK descriptor array for distributed matrix
# This array is identical to 
base.descinit <- function(dim, bldim, ldim, ICTXT=0)
{
  desc <- integer(9)
  desc[1L] <- 1L                    # matrix type
  desc[2L] <- ICTXT                 # CTXT_A
  desc[3L] <- dim[1L]               # M_A
  desc[4L] <- dim[2L]               # N_A
  desc[5L] <- bldim[1L]             # MB_A
  desc[6L] <- bldim[2L]             # NB_A
  desc[7L] <- 0L                    # RSRC_A
  desc[8L] <- 0L                    # CSRC_A
  desc[9L] <- max(1L, ldim[1L])     # LLD_A
  
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
    MYDIST <- (PROCS[i] + MYP[i] - ISRCPROC) %% PROCS[i]
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
      ldim[which(is.na(ldim))] <- 0L
    if (any(ldim<1)) ldim <- c(1L, 1L) # FIXME
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


# Hook into ScaLAPACK tool PDLAPRNT
base.rpdlaprnt <- function(m, n, a, desca)
{
  if (!is.double(a))
    storage.mode(a) <- "double"
  
  .Call("R_PDLAPRNT", 
        as.integer(m), as.integer(n),
        a, as.integer(desca),
        as.character(deparse(substitute(a))),
        6L,  #WCC: 0 for stderr, 6 for stdout. Both are disabled.
        PACKAGE="pbdBASE"
        )
  
  return( invisible(0) )
}

# Compute maximum dimension across all nodes
base.maxdim <- function(dim)
{
  mdim <- numeric(2)
  mdim[1] <- pbdMPI::allreduce(dim[1], op='max')
  mdim[2] <- pbdMPI::allreduce(dim[2], op='max')
  
  return( mdim )
}

# Compute dimensions on process MYROW=MYCOL=0
base.dim0 <- function(dim, ICTXT=0)
{
  blacs_ <- base.blacs(ICTXT=ICTXT)
  MYROW <- blacs_$MYROW
  MYCOL <- blacs_$MYCOL
  
  if (MYROW == 0 && MYCOL == 0){
    mx01 <- dim[1]
    mx02 <- dim[2]
  }
  
  mx01 <- pbdMPI::bcast(mx01)
  mx02 <- pbdMPI::bcast(mx02)
  
#  pbdMPI::barrier()
  
  if (MYROW==0 && MYCOL==0)
    return( dim )
  else
    return( c(mx01, mx02) )
}



# l2g and g2l
base.g2l_coord <- function(ind, dim, bldim, ICTXT=0)
{
  blacs_ <- base.blacs(ICTXT=ICTXT)
  procs <- c(blacs_$NPROW, blacs_$NPCOL)
  src <- c(0,0)
  
  out <- .Call("g2l_coords", 
                ind=as.integer(ind), dim=as.integer(dim), bldim=as.integer(bldim),
                procs=as.integer(procs), src=as.integer(src),
                PACKAGE="pbdBASE")
  
#  out[5:6] <- out[5:6] + 1
  
  if (out[3]!=blacs_$MYROW || out[4]!=blacs_$MYCOL)
    out <- rep(NA, 6)
  
  # out is a 'triple of pairs' stored as a length-6 vector, consisting of:
    # block position
    # process grid block
    # local coordinates
  # out will be a length 6 vector of NA when that global coord is not
  # relevant to the local storage
  
  return(out)
}

g2l_coord <- base.g2l_coord


base.l2g_coord <- function(ind, dim, bldim, ICTXT=0)
{
  blacs_ <- base.blacs(ICTXT=ICTXT)
  procs <- c(blacs_$NPROW, blacs_$NPCOL)
  myproc <- c(blacs_$MYROW, blacs_$MYCOL)
  
  out <- .Call("l2g_coords", 
                ind=as.integer(ind), dim=as.integer(dim), bldim=as.integer(bldim),
                procs=as.integer(procs), src=as.integer(myproc),
                PACKAGE="pbdBASE")
  
  return(out)
}

l2g_coord <- base.l2g_coord
