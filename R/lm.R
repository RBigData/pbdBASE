# Work in progress --- not ready for the big time

# lm.fit()
# Wrapper for custom PDGELS function, which solves linear least
# squares problems.
base.rpdgels <- function(a, b, tol=1e-7)
{
  # Matrix descriptors
  desca <- base.descinit(dim=a@dim, bldim=a@bldim, ldim=a@ldim, ICTXT=a@CTXT)
  descb <- base.descinit(dim=b@dim, bldim=b@bldim, ldim=b@ldim, ICTXT=b@CTXT)

  m <- desca[3]
  n <- desca[4]
  nrhs <- descb[4]
  
  # FIXME adjustment for weird lda issue
  mxlda <- pbdMPI::allreduce(desca[9], op='max')
  mxldb <- pbdMPI::allreduce(descb[9], op='max')
  
  desca[9] <- mxlda
  descb[9] <- mxldb
  
  IJ <- 1 # IA, JA, IB, JB
  
  # Determine size of work array 
  lwork <- .Fortran("RPDGELS", 
                    as.double(tol), as.character("N"), 
                    as.integer(m), as.integer(n), as.integer(nrhs),
                    double(1), as.integer(1), as.integer(1), as.integer(desca),
                    double(1), as.integer(1), as.integer(1), as.integer(descb),
                    double(1), double(1),
                    TAU=double(1), WORK=double(1), as.integer(-1), as.integer(1), 
                    integer(1), INFO=as.integer(0)
                    )$WORK[1]

  # Convert to .Call()
  # in: tol, trans, m, n, nrhs, IA, JA, desca, IB, JB, descb, lwork
  # in/out: A, B, 
  # out: FT, RSD, TAU, IPIV, RANK
  # local (just allocate in C, not a SEXP): work
  
  # IMPORTANT VVVVV
  # NOTE that FT must be initialized to zero
  # IMPORTANT ^^^^^

  # Fit the model
  out <- .Fortran("RPDGELS", 
                  TOL=as.double(tol), TRANS=as.character("N"),
                  M=as.integer(m), N=as.integer(n), NRHS=as.integer(nrhs),
                  A=a@Data, IA=as.integer(IJ), JA=as.integer(IJ), DESCA=as.integer(desca), 
                  B=b@Data, IB=as.integer(IJ), JB=as.integer(IJ), DESCB=as.integer(descb),
                  FT=matrix(double(prod(b@ldim)), b@ldim[1]), RSD=matrix(double(prod(b@ldim)), b@ldim[1]),
                  TAU=double(min(m, n)), WORK=double(lwork), LWORK=as.integer(lwork), IPIV=integer(a@ldim[2]),
                  RANK=integer(1), INFO=as.integer(0),
                  PACKAGE="pbdBASE")

  if (out$INFO!=0)
    warning(paste("ScaLAPACK returned INFO=", out$INFO, "; returned solution is likely invalid", sep=""))

  a@Data <- out$A
  b@Data <- out$B
  
  fitted.values <- new("ddmatrix", Data=out$FT, dim=b@dim,
                       ldim=b@ldim, bldim=b@bldim, CTXT=b@CTXT)
  residuals <- new("ddmatrix", Data=out$RSD, dim=b@dim,
                       ldim=b@ldim, bldim=b@bldim, CTXT=b@CTXT)
  
  # rearranging solution in the overdetermined and/or rank deficient case
  temp <- 1:n # indexing of coefficients
  b <- b[temp, ]
  
  
  # convert IPIV to global vector if it isn't already
  if (base.blacs(ICTXT=a@CTXT)$NPCOL > 1){
    c <- new("ddmatrix", Data=matrix(out$IPIV, nrow=1),
              dim=c(1, b@dim[1]), ldim=c(1, b@ldim[1]), 
              bldim=b@bldim, CTXT=b@CTXT)
    out$IPIV <- as.vector(c)
  }
  
  if (out$RANK < n){
    vec <- as.ddmatrix(matrix(NA, nrow=1, ncol=nrhs), bldim=b@bldim)
    b[(out$RANK+1):n, ] <- NA
    if (any(out$IPIV - temp != 0)){
      perm <- sapply(temp, function(i) temp[which(i==out$IPIV)])
      b <- b[perm, ]
    }
  }
  
  # rownames
#  if (base.ownany(dim=b@dim, bldim=b@bldim, CTXT=b@CTXT)){
#    coords <- sapply(temp, function(i) base.g2l_coord(ind=i, dim=b@dim, bldim=b@bldim, ICTXT=b@CTXT)[5])
#    mycoords <- coords[which(!is.na(coords))]
#    
#    rownames(b@Data) <- paste("x", mycoords, sep="")
#  }
  
  qr <- list(qr=a, tau=out$TAU, pivot=out$IPIV, tol=tol, rank=out$RANK)
  
  ret <- list(coefficients=b, residuals=residuals, effects=NULL, 
              rank=out$RANK, fitted.values=fitted.values, assign=NULL,
              qr=qr, df.residual=(a@dim[1] - out$RANK))
  
  return( ret )
}







# For n>=p case only.  In n<p case, have to do LQ
# qr()
base.rpdgeqrf <- function(x, tol=1e-7)
{
  # Matrix descriptors
  desca <- base.descinit(x@dim, x@bldim, x@ldim, ICTXT=x@CTXT)

  m <- desca[3]
  n <- desca[4]
  
  # Determine size of work array 
  lwork <- .Fortran("RPDGEQPF",
            TOL=as.double(tol), M=as.integer(m), N=as.integer(n), 
            A=double(1), as.integer(1), as.integer(1), DESCA=as.integer(desca), 
            IPIV=as.integer(1), TAU=double(1), 
            WORK=double(1), LWORK=as.integer(-1), 
            RANK=as.integer(n), INFO=integer(1),
            package="pbdBASE")$WORK[1]

  # perform QR
  out <- .Fortran("RPDGEQPF",
            TOL=as.double(tol), M=as.integer(m), N=as.integer(n), 
            A=x@Data, as.integer(1), as.integer(1), DESCA=as.integer(desca), 
            IPIV=integer(x@ldim[2]), TAU=double(min(m, n)), 
            WORK=double(lwork), LWORK=as.integer(lwork), 
            RANK=as.integer(n), INFO=integer(1),
            package="pbdBASE")

  if (out$INFO!=0)
    warning(paste("ScaLAPACK returned INFO=", out$INFO, "; returned solution is likely invalid", sep=""))

  x@Data <- out$A
  
  comm.print(out$TAU)

  ret <- list(qr=x, rank=out$RANK, tau=out$TAU, pivot=out$IPIV)
  
  attr(ret, "class") <- "qr"

  return( ret )
}








# qr.Q()
# recover Q from base.rpdgeqrf
base.pdorgqr <- function(qr)
{
  x <- qr$qr

  # Matrix descriptors
  desca <- base.descinit(x@dim, x@bldim, x@ldim, ICTXT=x@CTXT)

  m <- desca[3]
  n <- desca[4]
  
  k <- qr$rank
  
  # Determine size of work array 
  lwork <- .Fortran("PDORGQR",
            M=as.integer(m), N=as.integer(n), as.integer(k),
            A=double(1), as.integer(1), as.integer(1), DESCA=as.integer(desca), 
            TAU=double(1), 
            WORK=double(1), LWORK=as.integer(-1), INFO=integer(1),
            package="pbdBASE")$WORK[1]

  # perform QR
  out <- .Fortran("PDORGQR",
            M=as.integer(m), N=as.integer(n), as.integer(k),
            A=x@Data, as.integer(1), as.integer(1), DESCA=as.integer(desca), 
            TAU=as.double(qr$tau), 
            WORK=double(lwork), LWORK=as.integer(lwork), INFO=integer(1),
            package="pbdBASE")

  if (out$INFO!=0)
    warning(paste("ScaLAPACK returned INFO=", out$INFO, "; returned solution is likely invalid", sep=""))

  x@Data <- out$A
  
  return( x )
}





# qr.R
base.qr.R <- function(qr, complete=FALSE)
{
  ret <- qr$qr

  ret@Data <- base.low2zero(ret@Data, dim=ret@dim, bldim=ret@bldim)
  if (!complete)
    ret <- ret[1:min(ret@dim), ]
  
  return(ret)
}


