# ################################################
# ------------------------------------------------
# PDGESV:  Solving Ax=b
# ------------------------------------------------
# ################################################

base.rpdgesv <- function(a, b)
{
  ICTXT <- a@CTXT
  
  # Matrix descriptors
  desca <- base.descinit(a@dim, a@bldim, a@ldim, ICTXT=ICTXT)
  descb <- base.descinit(b@dim, b@bldim, b@ldim, ICTXT=ICTXT)
  
  n <- desca[4L]
  nrhs <- descb[4L]
  # max of the local dimensions
  mxldims <- c(base.maxdim(a@ldim), base.maxdim(b@ldim))
  
  # Call ScaLAPACK
  out <- .Call("R_PDGESV",
               as.integer(n), as.integer(nrhs), as.integer(mxldims),
               a@Data, as.integer(a@ldim), as.integer(desca),
               b@Data, as.integer(b@ldim), as.integer(descb),
               PACKAGE="pbdBASE"
              )
  
  if (out$info!=0)
    warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))
  
  return(out$B) 
}

# ################################################
# ------------------------------------------------
# PDGETRI:  Matrix inverse
# ------------------------------------------------
# ################################################

base.rpdgetri <- function(a)
{
  desca <- base.descinit(a@dim, a@bldim, a@ldim, ICTXT=a@CTXT)
  
  n <- desca[4]
  
#  lwork <- a@ldim[2] * a@bldim[2]
  
#  if (NPROW==NPCOL)
#    liwork <- a@ldim[2] + a@bldim[2]
#  else
#    liwork <- a@ldim[2] + 
#      max(ceiling(ceiling(a@ldim[1]/a@bldim[1])/(2/NPROW)), a@bldim[2])
  
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
# PDGESVD:  SVD of x
# ------------------------------------------------
# ################################################

base.rpdgesvd <- function(A, nu, nv)
{
  ICTXT <- A@CTXT
  
  # Matrix descriptors
  m <- A@dim[1]
  n <- A@dim[2]
  size <- min(A@dim)
  bldim <- A@bldim
  desca <- base.descinit(A@dim, A@bldim, A@ldim, ICTXT=ICTXT)

  if (nu==0){
    jobu <- 'N'
    udim <- c(1, 1)
  }
  else {
    jobu <- 'V'
    udim <- c(m, size)
  }
  if (nv==0){
    jobvt <- 'N'
    vtdim <- c(1,1)
  }
  else {
    jobvt <- 'V'
    vtdim <- c(size, n)
  }

  uldim <- base.numroc(udim, bldim, ICTXT=ICTXT)

  u <- new("ddmatrix", Data=matrix(nrow=0, ncol=0),
                       dim=udim, ldim=uldim, bldim=bldim, CTXT=ICTXT)
  descu <- base.descinit(u@dim, u@bldim, u@ldim, ICTXT=ICTXT)
  
  vtldim <- base.numroc(vtdim, bldim, ICTXT=ICTXT)

  vt <- new("ddmatrix", Data=matrix(nrow=0, ncol=0),
                        dim=vtdim, ldim=vtldim, bldim=bldim, CTXT=ICTXT)
  descvt <- base.descinit(vt@dim, vt@bldim, vt@ldim, ICTXT=ICTXT)

  mxa <- pbdMPI::allreduce(max(A@ldim), op='max')
  mxu <- pbdMPI::allreduce(max(uldim), op='max')
  mxvt <- pbdMPI::allreduce(max(vtldim), op='max')

  if (all(A@ldim==1))
    desca[9] <- mxa
  if (all(uldim==1))
    descu[9] <- mxu
  if (all(vtldim==1))
    descvt[9] <- mxvt

  if (A@dim[1]>1){
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

  # Call ScaLAPACK
  out <- .Call("R_PDGESVD", 
            as.integer(m), as.integer(n), as.integer(size),
            A@Data, as.integer(desca), as.integer(A@ldim),
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

  if (nu && nu < u@dim[2])
    u <- u[, 1L:nu]
  if (nv && nv < vt@dim[1])
    vt <- vt[1L:nv, ]

  if (out$info!=0)
    warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))

  ret <- list( d=out$d, u=u, vt=vt )

  return( ret )
}

# ################################################
# ------------------------------------------------
# PDGETRF:  LU Decomposition
# ------------------------------------------------
# ################################################

base.rpdgetrf <- function(a)
{
  # BLACS stuff
  ICTXT <- a@CTXT
  blacs_ <- base.blacs(ICTXT=ICTXT)
  MYROW <- blacs_$MYROW
  MYCOL <- blacs_$MYCOL

  desca <- base.descinit(a@dim, a@bldim, a@ldim, ICTXT=ICTXT)

  m <- desca[3]
  n <- desca[4]

  lipiv <- base.maxdim(a@ldim)[1] + a@bldim[1]

  # Call ScaLAPACK
  out <- .Call("R_PDGETRF",
               a@Data, as.integer(dim(a@Data)),
               as.integer(ICTXT), as.integer(MYROW), as.integer(MYCOL),
               as.integer(desca), as.integer(m), as.integer(n),
               as.integer(lipiv),
               PACKAGE="pbdBASE"
              )
  
  if (out$info!=0)
    warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))

  a@Data <- out$A
  return(a) 
}

# ################################################
# ------------------------------------------------
# PDPOTRF:  Cholesky Factorization
# ------------------------------------------------
# ################################################

base.pdpotrf <- function(a)
{
  desca <- base.descinit(a@dim, a@bldim, a@ldim, ICTXT=a@CTXT)
    
  n <- desca[4]
  
  uplo <- "U"
  
  # Call ScaLAPACK
  out <- .Call("R_PDPOTRF",
               as.integer(n),
               a@Data, as.integer(a@ldim), as.integer(desca),
               as.character(uplo),
               PACKAGE="pbdBASE"
              )
  
  if (out$info!=0)
    warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))
  else
    ret <- base.low2zero(A=out$A, dim=a@dim, ldim=a@ldim, bldim=a@bldim, CTXT=a@CTXT)
  
  return(ret) 
}

