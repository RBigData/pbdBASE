# ################################################
# ------------------------------------------------
# PDGESV:  Solving Ax=b
# ------------------------------------------------
# ################################################

base.rpdgesv <- function(a, b)
{
  # BLACS stuff
  ICTXT <- a@CTXT
  blacs_ <- base.blacs(ICTXT=ICTXT)
  MYROW <- blacs_$MYROW
  MYCOL <- blacs_$MYCOL

  # Matrix descriptors
  desca <- base.descinit(a@dim, a@bldim, a@ldim, ICTXT=ICTXT)
  descb <- base.descinit(b@dim, b@bldim, b@ldim, ICTXT=ICTXT)

  n <- desca[4L]
  nrhs <- descb[4L]
  mxldims <- c(base.maxdim(a@ldim), base.maxdim(b@ldim))

  # Call ScaLAPACK
  out <- .Call("R_PDGESV",
               a@Data, as.integer(dim(a@Data)),
               b@Data, as.integer(dim(b@Data)),
               as.integer(ICTXT), as.integer(MYROW), as.integer(MYCOL),
               as.integer(desca), as.integer(descb),
               as.integer(n), as.integer(nrhs),
               as.integer(mxldims),
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
  # BLACS stuff
  ICTXT <- a@CTXT
  blacs_ <- base.blacs(ICTXT=ICTXT)
  MYROW <- blacs_$MYROW
  MYCOL <- blacs_$MYCOL

  desca <- base.descinit(a@dim, a@bldim, a@ldim, ICTXT=ICTXT)

  n <- desca[4]

  lengths <- .Call("R_PDGETRISZ",
                   as.integer(ICTXT), as.integer(MYROW), as.integer(MYCOL),
                   as.integer(desca), as.integer(n),
                   PACKAGE="pbdBASE"
                  )

  out <- .Call("R_PDGETRI",
               a@Data, as.integer(dim(a@Data)), as.integer(ICTXT),
               as.integer(MYROW), as.integer(MYCOL),
               as.integer(desca), as.integer(n),
               as.integer(lengths$TEMP), as.integer(lengths$ITEMP),
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
  # BLACS stuff
  ICTXT <- A@CTXT
  blacs_ <- base.blacs(ICTXT=ICTXT)
  MYROW <- blacs_$MYROW
  MYCOL <- blacs_$MYCOL

  # Matrix descriptors
  m <- A@dim[1]
  n <- A@dim[2]
  size <- min(A@dim)
  bldim <- A@bldim
  desca <- base.descinit(A@dim, A@bldim, A@ldim, ICTXT=ICTXT)

  if (nu==0){
    jobu <- 'N'
    m <- 1
  }
  else
    jobu <- 'V'
  if (nv==0){
    jobvt <- 'N'
    n <- 1
  }
  else
    jobvt <- 'V'

  # Returns:  d, u, vt
  d <- double(0)

  udim <- c(m, size)
  uldim <- base.numroc(udim, bldim, ICTXT=ICTXT)

  u <- new("ddmatrix", Data=matrix(nrow=0, ncol=0),
                       dim=udim, ldim=uldim, bldim=bldim, CTXT=ICTXT)
  descu <- base.descinit(u@dim, u@bldim, u@ldim, ICTXT=ICTXT)
  
  vtdim <- c(size, n)
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

  if (pbdMPI::allreduce(desca[9], op='max')==1 && A@dim[1]>1)
    desca[9] <- mxa
  if (pbdMPI::allreduce(descu[9], op='max')==1 && u@dim[1]>1)
    desca[9] <- mxu
  if (pbdMPI::allreduce(descvt[9], op='max')==1 && vt@dim[1]>1)
    desca[9] <- mxvt

  # Call ScaLAPACK
  LWORK <- .Call("R_PDGESVDSZ",
                 as.integer(A@dim[1]), as.integer(A@dim[2]), as.integer(size),
                 as.integer(ICTXT), as.integer(MYROW), as.integer(MYCOL),
                 as.integer(desca), as.integer(descu), as.integer(descvt),
                 as.character(jobu), as.character(jobvt),
                 PACKAGE="pbdBASE"
                )

  out <- .Call("R_PDGESVD",
               as.integer(A@dim[1]), as.integer(A@dim[2]), as.integer(size),
               as.integer(ICTXT), as.integer(MYROW), as.integer(MYCOL),
               A@Data, as.integer(desca), as.integer(A@ldim),
               as.integer(uldim), as.integer(descu),
               as.integer(vtldim), as.integer(descvt),
               as.integer(LWORK$TEMP), as.character(jobu), as.character(jobvt),
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

  return( list( d=out$d, u=u, vt=vt ) )
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
  # BLACS stuff
  ICTXT <- a@CTXT
  blacs_ <- base.blacs(ICTXT=ICTXT)
  MYROW <- blacs_$MYROW
  MYCOL <- blacs_$MYCOL

  desca <- base.descinit(a@dim, a@bldim, a@ldim, ICTXT=ICTXT)

  n <- desca[4]
  
  uplo <- "U"

  # Call ScaLAPACK
  out <- .Call("R_PDPOTRF",
               a@Data, as.integer(dim(a@Data)),
               as.integer(ICTXT), as.integer(MYROW), as.integer(MYCOL),
               as.integer(desca), as.integer(n),
               as.character(uplo),
               PACKAGE="pbdBASE"
              )
  

  if (out$info!=0)
    warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))
  else
      out$A <- base.low2zero(A=out$A, dim=a@dim, ldim=a@ldim, bldim=a@bldim, CTXT=a@CTXT)

  return(out$A) 
}

