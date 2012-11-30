# ################################################
# ------------------------------------------------
# PDGESV:  Solving Ax=b
# ------------------------------------------------
# ################################################

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

# ################################################
# ------------------------------------------------
# PDGETRI:  Matrix inverse
# ------------------------------------------------
# ################################################

base.rpdgetri <- function(a)
{
  desca <- base.descinit(dim=a@dim, bldim=a@bldim, ldim=a@ldim, ICTXT=a@CTXT)
  
  n <- desca[4L]
  
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

# ################################################
# ------------------------------------------------
# PDGETRF:  LU Decomposition
# ------------------------------------------------
# ################################################

base.rpdgetrf <- function(a)
{
  desca <- base.descinit(dim=a@dim, bldim=a@bldim, ldim=a@ldim, ICTXT=a@CTXT)

  m <- desca[3L]
  n <- desca[4L]

  lipiv <- base.maxdim(a@ldim)[1L] + a@bldim[1L]

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

# ################################################
# ------------------------------------------------
# PDPOTRF:  Cholesky Factorization
# ------------------------------------------------
# ################################################

base.rpdpotrf <- function(x)
{
  desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@CTXT)
    
  n <- desca[4L]
  
  uplo <- "U"
  
  # Call ScaLAPACK
  out <- .Call("R_PDPOTRF",
               as.integer(n),
               x@Data, as.integer(x@ldim), as.integer(desca),
               as.character(uplo),
               PACKAGE="pbdBASE"
              )
  
  if (out$info!=0)
    warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))
  
  ret <- base.low2zero(A=out$A, dim=x@dim, ldim=x@ldim, bldim=x@bldim, CTXT=x@CTXT)
  
  return(ret) 
}

