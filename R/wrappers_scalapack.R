# ################################################
# ------------------------------------------------
# Linear Equations
# ------------------------------------------------
# ################################################

# ------------------------------------------------
# PDGESV:  Solving Ax=b
# ------------------------------------------------

base.rpdgesv <- function(a, b)
{
  ICTXT <- a@CTXT
  
  # Matrix descriptors
  desca <- base.descinit(dim=a@dim, bldim=a@bldim, ldim=a@ldim, ICTXT=ICTXT)
  descb <- base.descinit(dim=b@dim, bldim=b@bldim, ldim=b@ldim, ICTXT=ICTXT)
  
  n <- desca[4L]
  nrhs <- descb[4L]
  # max of the local dimensions
  mxldims <- c(base.maxdim(a@ldim), base.maxdim(b@ldim))

  if (!is.double(a@Data))
    storage.mode(a@Data) <- "double"
  if (!is.double(b@Data))
    storage.mode(b@Data) <- "double"

  # Call ScaLAPACK
  out <- .Call("R_PDGESV",
               as.integer(n), as.integer(nrhs), as.integer(mxldims),
               a@Data, as.integer(a@ldim), as.integer(desca),
               b@Data, as.integer(b@ldim), as.integer(descb),
               PACKAGE="pbdBASE"
              )
  
  if (out$info!=0)
    warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))
  
  b@Data <- out$B
  
  return(b) 
}

# ------------------------------------------------
# PDGETRI:  Matrix inverse
# ------------------------------------------------

base.rpdgetri <- function(a)
{
  desca <- base.descinit(dim=a@dim, bldim=a@bldim, ldim=a@ldim, ICTXT=a@CTXT)
  
  n <- desca[4L]
  
  if (!is.double(a@Data))
    storage.mode(a@Data) <- "double"
  
  out <- .Call("R_PDGETRI",
               a@Data, as.integer(a@ldim), 
               as.integer(desca), as.integer(n),
               PACKAGE="pbdBASE"
              )
  
  if (out$info!=0)
    warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))

  return(out$A) 
}

# ################################################
# ------------------------------------------------
# Matrix Factorizations
# ------------------------------------------------
# ################################################

# ------------------------------------------------
# PDGESVD:  SVD of x
# ------------------------------------------------

base.rpdgesvd <- function(x, nu, nv)
{
  ICTXT <- x@CTXT
  
  # Matrix descriptors
  m <- x@dim[1L]
  n <- x@dim[2L]
  size <- min(x@dim)
  bldim <- x@bldim
  desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=ICTXT)

  if (nu==0){
    jobu <- 'N'
    udim <- c(1L, 1L)
  }
  else {
    jobu <- 'V'
    udim <- c(m, size)
  }
  if (nv==0){
    jobvt <- 'N'
    vtdim <- c(1L, 1L)
  }
  else {
    jobvt <- 'V'
    vtdim <- c(size, n)
  }

  uldim <- base.numroc(dim=udim, bldim=bldim, ICTXT=ICTXT)

  u <- new("ddmatrix", Data=matrix(nrow=0, ncol=0),
                       dim=udim, ldim=uldim, bldim=bldim, CTXT=ICTXT)
  descu <- base.descinit(dim=u@dim, bldim=u@bldim, ldim=u@ldim, ICTXT=ICTXT)
  
  vtldim <- base.numroc(dim=vtdim, bldim=bldim, ICTXT=ICTXT)

  vt <- new("ddmatrix", Data=matrix(nrow=0, ncol=0),
                        dim=vtdim, ldim=vtldim, bldim=bldim, CTXT=ICTXT)
  descvt <- base.descinit(dim=vt@dim, bldim=vt@bldim, ldim=vt@ldim, ICTXT=ICTXT)

  mxa <- pbdMPI::allreduce(max(x@ldim), op='max')
  mxu <- pbdMPI::allreduce(max(uldim), op='max')
  mxvt <- pbdMPI::allreduce(max(vtldim), op='max')

  if (all(x@ldim==1))
    desca[9] <- mxa
  if (all(uldim==1))
    descu[9] <- mxu
  if (all(vtldim==1))
    descvt[9] <- mxvt

  if (x@dim[1]>1){
    if (pbdMPI::allreduce(desca[9], op='max')==1)
      desca[9] <- mxa
  }
  if (u@dim[1]>1){
    if (pbdMPI::allreduce(descu[9], op='max')==1)
      desca[9] <- mxu
  }
  if (vt@dim[1]>1){
    if (pbdMPI::allreduce(descvt[9], op='max')==1)
      desca[9] <- mxvt
  }

  if (!is.double(x@Data))
    storage.mode(x@Data) <- "double"

  # Call ScaLAPACK
  out <- .Call("R_PDGESVD", 
            as.integer(m), as.integer(n), as.integer(size),
            x@Data, as.integer(desca), as.integer(x@ldim),
            as.integer(uldim), as.integer(descu),
            as.integer(vtldim), as.integer(descvt),
            as.character(jobu), as.character(jobvt),
            PACKAGE="pbdBASE"
          )

  if (nu==0)
    u <- NULL
  else 
    u@Data <- out$u
  if (nv==0)
    vt <- NULL
  else
    vt@Data <- out$vt

  if (nu && nu < u@dim[2L])
    u <- u[, 1L:nu]
  if (nv && nv < vt@dim[1L])
    vt <- vt[1L:nv, ]

  if (out$info!=0)
    warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))

  ret <- list( d=out$d, u=u, vt=vt )

  return( ret )
}

# ------------------------------------------------
# PDGETRF:  LU Decomposition
# ------------------------------------------------

base.rpdgetrf <- function(a)
{
  desca <- base.descinit(dim=a@dim, bldim=a@bldim, ldim=a@ldim, ICTXT=a@CTXT)

  m <- desca[3L]
  n <- desca[4L]

  lipiv <- base.maxdim(a@ldim)[1L] + a@bldim[1L]

  if (!is.double(a@Data))
    storage.mode(a@Data) <- "double"

  # Call ScaLAPACK
  out <- .Call("R_PDGETRF",
               as.integer(m), as.integer(n),
               a@Data, as.integer(dim(a@Data)), as.integer(desca),
               as.integer(lipiv),
               PACKAGE="pbdBASE"
              )
  
  if (out$info!=0)
    warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))

  a@Data <- out$A
  return(a) 
}

# ------------------------------------------------
# PDPOTRF:  Cholesky Factorization
# ------------------------------------------------

base.rpdpotrf <- function(x)
{
  desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@CTXT)
    
  n <- desca[4L]
  
  uplo <- "U"
  
  if (!is.double(x@Data))
    storage.mode(x@Data) <- "double"
  
  # Call ScaLAPACK
  out <- .Call("R_PDPOTRF",
               as.integer(n),
               x@Data, as.integer(x@ldim), as.integer(desca),
               as.character(uplo),
               PACKAGE="pbdBASE"
              )
  
  if (out$info!=0)
    warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))
  
  ret <- new("ddmatrix", Data=out$A, dim=x@dim, bldim=x@bldim, CTXT=x@CTXT)
  
#  ret <- base.tri2zero(dx=ret, 'L', 'N')
  
  return(ret) 
}



# ################################################
# ------------------------------------------------
# Auxillary
# ------------------------------------------------
# ################################################

base.indxg2p <- function(INDXGLOB, NB, NPROCS)
{
  
  ISRCPROC <- 0L
  
  ret <- (ISRCPROC + (INDXGLOB - 1L) / NB) %% NPROCS
  
  return( ret )
}


numroc2 <- function(N, NB, IPROC, NPROCS)
{
  ISRCPROC <- 0L
  
  MYDIST <- (NPROCS + IPROC -  ISRCPROC) %% NPROCS
  NBLOCKS <- floor(N / NB)
  ldim <- floor(NBLOCKS / NPROCS) * NB
  EXTRABLKS <- NBLOCKS %% NPROCS
  
  if (is.na(EXTRABLKS))
    EXTRABLKS <- 0L
  
  if (MYDIST < EXTRABLKS)
    ldim <- ldim + NB
  else if (MYDIST == EXTRABLKS)
    ldim <- ldim + N %% NB
  
  return(ldim)
}


# matrix norms
base.rpdlange <- function(x, type)
{
  desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@CTXT)
  
  m <- x@dim[1L]
  n <- x@dim[2L]
  
  if (length(type)>1L)
    type <- type[1L]
  
  type <- toupper(type)
  
#  if (type == "M" || type == "F")
#    lwork <- 1L
#  else {
#    blacs_ <- base.blacs(ICTXT=x@CTXT)
#    ia <- ja <- 1L
#    
#    if (type == "O"){
#      mb <- x@bldim[1L]
#      npcol <- blacs_$NPCOL
#      
#      icoffa <- (ja-1L) %% mb
#      iacol <- base.indxg2p(ja, mb, npcol)
#      
#      lwork <- numroc2(n+icoffa, mb, blacs_$MYCOL, npcol)
#    }
#    else if (type == "I"){
#      nb <- x@bldim[2L]
#      nprow <- blacs_$NPROW
#      
#      iroffa <- (ia-1L) %% nb
#      iarow <- base.indxg2p(ia, nb, nprow)
#      
#      lwork <- numroc2(m+iroffa, nb, blacs_$MYROW, nprow)
#    }
#  }
#  
#  lwork <- max(lwork, 1L)
  
  if (!is.double(x@Data))
    storage.mode(x@Data) <- "double"
  
  ret <- .Call("R_PDLANGE", 
        as.character(type), as.integer(m), as.integer(n),
        x@Data, as.integer(desca),
        PACKAGE="pbdBASE")
  
  return( ret )
}


# Inverse condition number - general matrix
base.rpdgecon <- function(x, type)
{
  desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@CTXT)
  
  m <- x@dim[1L]
  n <- x@dim[2L]
  
  if (length(type)>1L)
    type <- type[1L]
  
  type <- toupper(type)
  
  if (!is.double(x@Data))
    storage.mode(x@Data) <- "double"
  
  ret <- .Call("R_PDGECON", 
        as.character(type), as.integer(m), as.integer(n),
        x@Data, as.integer(desca), as.integer(x@ldim),
        PACKAGE="pbdBASE")
  
  if (ret[2] < 0)
    warning(paste("INFO =", ret[2]))
  
  return( ret[1] )
}


# Inverse condition number - triangular matrix
base.rpdtrcon <- function(x, type, uplo="L")
{
  desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@CTXT)
  
#  m <- x@dim[1L]
  n <- x@dim[2L]
  
  if (length(type)>1L)
    type <- type[1L]
  
  type <- toupper(type)
  uplo <- toupper(uplo)
  
  if (!is.double(x@Data))
    storage.mode(x@Data) <- "double"
  
  ret <- .Call("R_PDTRCON", 
        as.character(type), as.character(uplo), 
        as.integer(n), x@Data, as.integer(desca),
        PACKAGE="pbdBASE")
  
  if (ret[2] < 0)
    warning(paste("INFO =", ret[2]))
  
  return( ret[1] )
}




