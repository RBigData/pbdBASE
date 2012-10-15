# dropper function, used in subsetting
base.dropper <- function(x, oldbldim, iorj, ij, ICTXT)
{
  blacs_ <- base.blacs(ICTXT)

  bldim <- oldbldim #x@bldim
  if (x@CTXT != ICTXT){
# FIXME: would like to alter block dimension so that the data
# rebalances, but the code below causes a horrible null pointer
# problem for some cases.
#    if (ICTXT==1)
#      bldim <- c(dim(x)[1], ceiling(bldim[2] / blacs_$NPCOL))
#    if (ICTXT==2)
#      bldim <- c(ceiling(bldim[1] / blacs_$NPROW), dim(x)[2])

    newObj <- base.reblock(dx=x, bldim=bldim, ICTXT)
  }

#  if (blacs_$MYROW != -1 && blacs_$MYCOL != -1){
    if (iorj=='i'){ # rows
      if (newObj@ldim[1] == newObj@dim[1]){
        new <- newObj@Data[ij, ]
        if (base::length(new)==0)
          new <- matrix(0)
        if (!is.matrix(new))
          new <- matrix(new, ncol=newObj@ldim[2])

        newObj@Data <- new
      }
    } else { # columns
      if (newObj@ldim[2] == newObj@dim[2]){
        new <- newObj@Data[, ij]
        if (base::length(new)==0)
          new <- matrix(0)
        if (!is.matrix(new))
          new <- matrix(new, nrow=newObj@ldim[1])

        newObj@Data <- new
      }
    }
#  }

  if (iorj=='i'){
    if (base::any(ij<0))
      newObj@dim[1] <- newObj@dim[1] - base::length(ij)
    else
      newObj@dim[1] <- base::length(ij)
  } else {
    if (base::any(ij<0))
      newObj@dim[2] <- newObj@dim[2] - base::length(ij)
    else
      newObj@dim[2] <- base::length(ij)
  }

  newObj@ldim <- dim(newObj@Data)
  
  return(newObj)
}

dropper <- base.dropper


# For use with local arithmetic; basically does nothing if the 
# local storage is just filler to make scalapack happy
# return is logical answer to the question:  'do I own anything?', 
# and not a C-style return
base.ownany <- function(dim, bldim, CTXT=0)
{
  if (length(bldim)==1)
    bldim <- rep(bldim, 2)
  check <- base.numroc(dim=dim, bldim=bldim, ICTXT=CTXT, fixme=FALSE)
  
  if (any(check<1))
    return(FALSE)
  else
    return(TRUE)
}

ownany <- function(x, ..., dim, bldim, CTXT=0)
{
  if (!missing(bldim) && length(bldim)==1)
    bldim <- rep(bldim, 2)
  if (!missing(x) && is.ddmatrix(x))
    return( base.ownany(dim=x@dim, bldim=x@bldim, CTXT=x@CTXT) )
  else if (!missing(dim) && !missing(bldim) && missing(x) && is.numeric(dim) && is.numeric(bldim))
    return( base.ownany(dim=dim, bldim=bldim, CTXT=CTXT) )
  else{
    print("Error: bad input(s) in ownany()")
    stop("")
  }
}


# checking compatibility between distributed matrices for use with
# scalapack/pblas. For internal use only.
base.checkem <- function(x, y, checks=1:3)
{
  # All dimension equal
  if (1 %in% checks)
    if (any(x@dim!=y@dim)){
      pbdMPI::comm.print("Error: non-conformable distributed arrays")
      stop("")
    }
  # Same BLACS context
  if (2 %in% checks)
    if (x@CTXT != y@CTXT){
      pbdMPI::comm.print("Error: Distributed matrices 'x' and 'y' must belong to the same BLACS context")
      stop("")
    }
  # Same blocking dimension
  if (3 %in% checks)
    if (any(x@bldim != y@bldim)){
      pbdMPI::comm.print("Distributed matrices 'x' and 'y' must have the same block dimension.")
      stop("")
    }
}

checkem <- base.checkem


# Compute maximum dimension across all nodes
base.maxdim <- function(dim)
{
  mdim <- numeric(2)
  mdim[1] <- pbdMPI::allreduce(dim[1], op='max')
  mdim[2] <- pbdMPI::allreduce(dim[2], op='max')
  
  return(mdim)
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



# Take a (regular) matrix common to all nodes and distribute it as
# ddmatrix.  Should only be used in testing! This is inefficient 
# for real work.
base.submat <- function(A, bldim, ICTXT=0)
{
  blacs_ <- base.blacs(ICTXT=ICTXT)
  myP <- c(blacs_$MYROW, blacs_$MYCOL)
  PROCS <- c(blacs_$NPROW, blacs_$NPCOL)
  RSRC <- CSRC <- 0 # processes with first row/col of global A
  ISRCPROC <- 0
  
  if (length(bldim)==1) bldim <- rep(bldim, 2)
  dim <- dim(A)
  ldim <- base.numroc(dim, bldim)
  
  out <- .Call("block_submat", 
             A=A,
             dim=as.integer(dim),
             ldim=as.integer(ldim),
             bldim=as.integer(bldim),
             gP=as.integer(PROCS),
             myP=as.integer(myP),
             SRC=as.integer(c(RSRC, CSRC)),
             PACKAGE="pbdBASE"
             )
  
  return(out)
}

# Reverse of submat above.  Same restrictions apply.
base.gmat <- function(dx, proc.dest="all")
{
  blacs_ <- base.blacs(dx@CTXT)
  myP <- c(blacs_$MYROW, blacs_$MYCOL)
  PROCS <- c(blacs_$NPROW, blacs_$NPCOL)
  RSRC <- CSRC <- 0 # processes with first row/col of global A
  ISRCPROC <- 0
  
  xattrs <- attributes(dx@Data)
  
  dim <- dx@dim
  ldim <- dx@ldim
  bldim <- dx@bldim
  
  if (any(dim==0)){
    if (proc.dest=="all" || proc.dest==comm.rank())
      out <- matrix(nrow=dim[1], ncol=dim[2])
    else
      out <- NULL
    return(out)
  }
  else  
    out <- .Call("submat_to_gmat", 
               subx=dx@Data,
               dim=as.integer(dim),
               bldim=as.integer(bldim),
               gP=as.integer(PROCS),
               myP=as.integer(myP),
               SRC=as.integer(c(RSRC, CSRC)),
               PACKAGE="pbdBASE"
               )

  if (all(proc.dest=="all"))
    out <- pbdMPI::allreduce(out, op='sum')
  else {
    if (all(myP==proc.dest))
      outproc <- pbdMPI::comm.rank()
    else
      outproc <- 0
    outproc <- pbdMPI::allreduce(outproc, op='max')
    out <- pbdMPI::reduce(out, op='sum', rank.dest=outproc)
  }
  
  if (is.null(out))
    return(out)
  else {
    out <- matrix(out, nrow=dim[1], ncol=dim[2])
    if (length(xattrs)>1){
      oattrs <- union(attributes(out), xattrs[-1])
      names(oattrs) <- names(xattrs)
      attributes(out) <- oattrs
    }
    return( out )
  }
}


# print first few entries of the global matrix
base.firstfew <- function(dx, atmost=5)
{
  blacs_ <- base.blacs(dx@CTXT)
  MYROW <- blacs_$MYROW
  MYCOL <- blacs_$MYCOL
  NPROW <- blacs_$NPROW
  NPCOL <- blacs_$NPCOL

  MB <- dx@bldim[1]
  NB <- dx@bldim[2]

  if (prod(dx@dim) < atmost)
    atmost <- prod(dx@dim)

  dim <- c( min(dx@dim[1], atmost), min(dx@dim[2], atmost) )

  out <- numeric(atmost)
  ct <- 1
  for (j in 1:dim[2]-1){
    for (i in 1:dim[1]-1){
      l <- floor(i / (NPROW * MB))
      m <- floor(j / (NPCOL * NB))
      
      pr <- (0 + floor(i/MB)) %% NPROW
      pc <- (0 + floor(j/NB)) %% NPCOL
      
      if (MYROW==pr && MYCOL==pc){
        x <- 1 + i %% MB;
        y <- 1 + j %% NB;
        out[ct] <- dx@Data[x+MB*l,y+NB*m]

        ct <- ct+1
      }
      ct <- pbdMPI::allreduce(ct, op='max')
      if (ct == atmost+1)
         break
      barrier()
    }
    if (ct == atmost+1)
      break
  }
  barrier()
  out <- pbdMPI::allreduce(out, op='sum')
  return(out)
}

# fill lower triangle of distributed matrix A with 0's
# used to force chol() return from ScaLAPACK to match R's return
base.low2zero <- function(A, dim, ldim, bldim, CTXT=0)
{
  blacs_ <- base.blacs(CTXT)
  myproc <- c(blacs_$MYROW, blacs_$MYCOL)
  PROCS <- c(blacs_$NPROW, blacs_$NPCOL)
  RSRC <- CSRC <- 0
  ISRCPROC <- 0
  
  out <- .Call("fill_lower_zero", 
             A=A,
             dim=as.integer(dim),
             bldim=as.integer(bldim),
             gP=as.integer(PROCS),
             myP=as.integer(myproc),
             SRC=as.integer(c(RSRC, CSRC)),
             PACKAGE="pbdBASE"
             )
  
  return(out) # A@Data
}

#base.reblock <- function(dx, bldim=dx@bldim, ICTXT)
#{
#  if (length(bldim)==1)
#    bldim <- rep(bldim, 2)

#  dim <- dx@dim
#  
#  xattrs <- attributes(dx@Data)
#  
#  m <- dim[1]
#  n <- dim[2]
#  
#  descx <- base.descinit(dim, dx@bldim, dx@ldim, ICTXT=dx@CTXT)
#  
#  # create descriptor for and initialize reblocked output matrix
#  ldimB <- base.numroc(dim, bldim, ICTXT=ICTXT)

#  # Otherwise
#  descb <- base.descinit(dim, bldim, ldimB, ICTXT=ICTXT)
#  dB <- new("ddmatrix", Data=matrix(0, ldimB[1], ldimB[2]), 
#           dim=dim, ldim=ldimB, bldim=bldim, CTXT=ICTXT)

#  xblacs_ <- base.blacs(dx@CTXT)
#  if (xblacs_$MYROW==-1 || xblacs_$MYCOL==-1){
##    descx <- rep(0, 9)
#    descx[2] <- -1
#  }

#  blacs_ <- base.blacs(ICTXT=ICTXT)
#  if (blacs_$MYROW==-1 || blacs_$MYCOL==-1){
##    descb <- rep(0, 9)
#    descb[2] <- -1
#  }
#  
#  # lda's of 1 infuriate pdgemr2d
#  mxx <- pbdMPI::allreduce(max(dx@ldim), op='max')
#  mxb <- pbdMPI::allreduce(max(ldimB), op='max')

#  if (all(ldim(dx)==1))
#    descx[9] <- mxx
#  if (all(ldim(dB)==1))
#    descb[9] <- mxb

##   things break when a 1xm matrix is stored so we have to pass a larger lda
#  if (pbdMPI::allreduce(as.integer(descx[9]), op='max')==1 && dx@dim[1]>1)
#    descx[9] <- mxx
#  if (pbdMPI::allreduce(as.integer(descb[9]), op='max')==1)
#    descb[9] <- mxb

#  ### WCC: dim(b@Data) is added to allocate B for replacing.
#  dB@Data <- .Call("R_PDGEMR2D",
#                   as.integer(m), as.integer(n),
#                   dx@Data, as.integer(descx),
#                   as.integer(ldimB), as.integer(descb),
#                   as.integer(0),
#                   as.integer(ldim(dx)[1]), as.integer(ldim(dx)[2]),
#                   as.integer(ldim(dB)[1]), as.integer(ldim(dB)[2]),
#                   PACKAGE="pbdBASE"
#                   )

#  if (length(xattrs)>1){
#    battrs <- union(attributes(dB@Data), xattrs[-1])
#    names(battrs) <- names(xattrs)
#    attributes(dB@Data) <- battrs
#  }

#  return(dB)
#}

base.reblock <- function(dx, bldim=dx@bldim, ICTXT)
{
  if (length(bldim)==1)
    bldim <- rep(bldim, 2)

#  blacs_ <- base.blacs(ICTXT=ICTXT)
#  ICTXT <- blacs_$ICTXT + 0
  
  dim <- dx@dim
  m <- dim[1]
  n <- dim[2]
  xattrs <- attributes(dx@Data)

  ldimB <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
  TldimB <- ldimB # true ldimB

  # lda's of 1 infuriate pdgemr2d
  mxx <- pbdMPI::allreduce(max(dx@ldim), op='max')
  mxb <- pbdMPI::allreduce(max(ldimB), op='max')

  if (all(dx@ldim==1))
    dx@ldim[1] <- mxx
  if (all(ldimB==1))
    ldimB[1] <- mxb

  if (pbdMPI::allreduce(dx@ldim[1], op='max')==1 && dx@dim[1]>1)
    dx@ldim[1] <- mxx
  if (pbdMPI::allreduce(ldimB[1], op='max')==1)
    ldimB[1] <- mxb

  descx <- base.descinit(dim=dim, bldim=dx@bldim, ldim=dx@ldim, ICTXT=dx@CTXT)
  descb <- base.descinit(dim=dim, bldim=bldim, ldim=ldimB, ICTXT=ICTXT)

  dB <- new("ddmatrix", Data=matrix(0, 1, 1), 
           dim=dim, ldim=TldimB, bldim=bldim, CTXT=ICTXT)

  xblacs_ <- base.blacs(dx@CTXT)
  if (xblacs_$MYROW==-1 || xblacs_$MYCOL==-1){
#    descx <- rep(0, 9)
    descx[2] <- -1
  }

  blacs_ <- base.blacs(ICTXT=ICTXT)
  if (blacs_$MYROW==-1 || blacs_$MYCOL==-1){
#    descb <- rep(0, 9)
    descb[2] <- -1
  }

  ret <- .Call("R_PDGEMR2D",
               as.integer(m), as.integer(n),
               dx@Data, as.integer(descx),
               as.integer(TldimB), as.integer(descb),
               as.integer(0), # context 0 is always passed since pdgemr2d 
               # requires the grids to have at least 1 processor in common
               as.integer(dx@ldim[1]), as.integer(dx@ldim[2]),
               PACKAGE="pbdBASE"
            )

#  ret <- .Fortran("PDGEMR2D",
#                   as.integer(m), as.integer(n),
#                   dx@Data, 
#                   as.integer(1), as.integer(1), as.integer(descx),
#                   B=matrix(0, TldimB[1], TldimB[2]),
#                   as.integer(1), as.integer(1), as.integer(descb),
#                   as.integer(ICTXT),
#                   PACKAGE="pbdBASE", DUP=T
#                  )$B

    ret <- ret + 0
    dB@Data <- ret
    
  if (length(xattrs) > 1){
    battrs <- union(attributes(dB@Data), xattrs[-1])
    names(battrs) <- names(xattrs)
    attributes(dB@Data) <- battrs
  }

  return(dB)
}

reblock <- base.reblock


# l2g and g2l
base.g2l_coord <- function(ind, dim, bldim, ICTXT=0)
{
  blacs_ <- base.blacs(ICTXT=ICTXT)
  procs <- c(blacs_$NPROW, blacs_$NPCOL)
  src <- c(0,0)
  
  out <- .Call("rcpp_g2l_coord", 
                ind=as.integer(ind),
                dim=as.integer(dim), bldim=as.integer(bldim),
                procs=as.integer(procs), src=as.integer(src),
                PACKAGE="pbdBASE"
               )
  
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
  myproc <- pbdMPI::comm.rank()
  
  out <- .Call("rcpp_l2g_coord", 
                ind=as.integer(ind),
                dim=as.integer(dim), bldim=as.integer(bldim),
                procs=as.integer(procs), src=as.integer(myproc),
                PACKAGE="pbdBASE"
               )
  
  return(out)
}

l2g_coord <- base.l2g_coord


base.mat.to.ddmat <- function(x, bldim=.BLDIM, ICTXT=0)
{
  if (!is.matrix(x)) {
    comm.print("input 'x' must be a matrix") 
    stop("")
  }
  else if (length(bldim) == 1) 
    bldim <- rep(bldim, 2) 
  else if (diff(bldim) != 0)
    warning("Most ScaLAPACK routines do not allow for non-square blocking.  This is highly non-advised.")

  blacs_ <- blacs(ICTXT=ICTXT)
  nprows <- blacs_$NPROW
  npcols <- blacs_$NPCOL
  dim <- dim(x)
  ldim <- base.numroc(dim, bldim)

  ddmat <- new("ddmatrix", Data=matrix(0), dim=dim, ldim=ldim, bldim=bldim)

  if (any(ldim==0)){
    ddmat@ldim <- c(1,1)
    ddmat@Data <- matrix(0)
  } else
      ddmat@Data <- base.submat(A=x, bldim=bldim)
  
#  pbdMPI::barrier()
  return(ddmat)
}

# wrapper around some C++ code to handle R's cyclic matrix-vector operatoins
# for distributed matrices, e.g. preserving
# matrix(1:6, ncol=2) + 1:2
base.vecops <- function(dx, vec, FUN)
{
  blacs_ <- base.blacs(dx@CTXT)
  procs <- c(blacs_$NPROW, blacs_$NPCOL)
  myprocs <- c(blacs_$MYROW, blacs_$MYCOL)
  src <- c(0,0)
  
  out <- .Call("ddmatrix_vecops", 
    dx@Data, as.integer(dx@dim), as.integer(dx@bldim),
    as.double(vec), as.integer(length(vec)),
    as.integer(procs), as.integer(myprocs), as.integer(src), 
    as.integer(FUN),
    PACKAGE="pbdBASE"
  )
  
  return(out)
}

# same as above, but just for insertion
base.insert <- function(dx, vec, i, j)
{
  blacs_ <- base.blacs(dx@CTXT)
  procs <- c(blacs_$NPROW, blacs_$NPCOL)
  myprocs <- c(blacs_$MYROW, blacs_$MYCOL)
  src <- c(0,0)
  
  if (all(i<0)){
    new <- 1:dx@dim[1]
    i <- new[-which(new %in% abs(i))] # FIXME make this less stupid
  }
  if (all(j<0)){
    new <- 1:dx@dim[2]
    j <- new[-which(new %in% abs(j))] # FIXME make this less stupid
  }
  
  out <- .Call("ddmatrix_insert", 
    dx@Data, as.integer(dx@dim), as.integer(dx@bldim),
    as.double(vec), as.integer(length(vec)),
    as.integer(i), as.integer(length(i)), as.integer(j), as.integer(length(j)),
    as.integer(procs), as.integer(myprocs), as.integer(src),
    PACKAGE="pbdBASE"
  )
  
  return(out)
}

