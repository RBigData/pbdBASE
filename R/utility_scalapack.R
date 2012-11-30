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


# Hook into ScaLAPACK tool PDLAPRNT
base.rpdlaprnt <- function(dx)
{
  m <- dx@dim[1L]
  n <- dx@dim[2L]
  
  desca <- base.descinit(dim=dx@dim, bldim=dx@bldim, ldim=dx@ldim, ICTXT=dx@CTXT)
  
  .Call("R_PDLAPRNT", 
        as.integer(m), as.integer(n),
        as.double(dx@Data), as.integer(desca),
        as.character(deparse(substitute(dx))),
        as.integer(6),  #WCC: 0 for stderr, 6 for stdout. Both are disabled.
        PACKAGE="pbdBASE"
        )
  
  return( invisible(0) )
}

rpdlaprnt <- base.rpdlaprnt
