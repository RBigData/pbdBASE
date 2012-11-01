# lm.fit()
# Wrapper for custom PDGELS function, which solves linear least
# squares problems.
base.rpdgels <- function(a, b, tol=1e-7)
{
  oldctxt <- b@CTXT
  
  # Matrix descriptors
  desca <- base.descinit(dim=a@dim, bldim=a@bldim, ldim=a@ldim, ICTXT=oldctxt)
  descb <- base.descinit(dim=b@dim, bldim=b@bldim, ldim=b@ldim, ICTXT=oldctxt)
  
  m <- desca[3]
  n <- desca[4]
  nrhs <- descb[4]
  
  # FIXME adjustment for weird lda issue
  mxlda <- pbdMPI::allreduce(desca[9], op='max')
  mxldb <- pbdMPI::allreduce(descb[9], op='max')
  
  if (desca[9]==1)
    desca[9] <- mxlda
  if (descb[9]==1)
    descb[9] <- mxldb
  
  IJ <- 1 # IA, JA, IB, JB
  
#  # Determine size of work array 
#  lwork <- .Fortran("RPDGELS", 
#                    as.double(tol), as.character("N"), 
#                    as.integer(m), as.integer(n), as.integer(nrhs),
#                    double(1), as.integer(IJ), as.integer(IJ), as.integer(desca),
#                    double(1), as.integer(IJ), as.integer(IJ), as.integer(descb),
#                    double(1), double(1),
#                    TAU=double(1), WORK=double(1), as.integer(-1), as.integer(1), 
#                    integer(1), INFO=as.integer(0)
#                    )$WORK[1]

#  # Convert to .Call()
#  # in: tol, trans, m, n, nrhs, IA, JA, desca, IB, JB, descb, lwork
#  # in/out: A, B, 
#  # out: FT, RSD, TAU, IPIV, RANK
#  # local (just allocate in C, not a SEXP): work
#  
#  # IMPORTANT VVVVV
#  # NOTE that FT must be initialized to zero
#  # IMPORTANT ^^^^^

#  # Fit the model
#  out <- .Fortran("RPDGELS", 
#                  TOL=as.double(tol), TRANS=as.character("N"),
#                  M=as.integer(m), N=as.integer(n), NRHS=as.integer(nrhs),
#                  A=a@Data, IA=as.integer(IJ), JA=as.integer(IJ), DESCA=as.integer(desca), 
#                  B=b@Data, IB=as.integer(IJ), JB=as.integer(IJ), DESCB=as.integer(descb),
#                  FT=matrix(double(prod(b@ldim)), b@ldim[1]), RSD=matrix(double(prod(b@ldim)), b@ldim[1]),
#                  TAU=double(min(m, n)), WORK=double(lwork), LWORK=as.integer(lwork), IPIV=integer(a@ldim[2]),
#                  RANK=integer(1), INFO=as.integer(0),
#                  PACKAGE="pbdBASE")

  out <- .Call("R_PDGELS",
            TOL=as.double(tol), M=as.integer(m), N=as.integer(n), NRHS=as.integer(nrhs),
            A=a@Data, ALDIM=as.integer(a@ldim), DESCA=as.integer(desca),
            B=b@Data, BLDIM=as.integer(b@ldim), DESCB=as.integer(descb),
            LTAU=as.integer(min(m, n)),
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
  temp <- 1L:n # indexing of coefficients
  if (m >= n){
    b <- b[temp, , ICTXT=1L]
  } else {
    cdim <- c(n-b@dim[1L], b@dim[2L])
    cldim <- base.numroc(dim=cdim, bldim=b@bldim, ICTXT=b@CTXT, fixme=TRUE)
    c <- new("ddmatrix", Data=matrix(as.double(NA), nrow=cldim[1L], ncol=cldim[2L]), dim=cdim, ldim=cldim, bldim=b@bldim, CTXT=b@CTXT)
    b <- base.rbind(b, c, ICTXT=1L)
  }
  
  # convert IPIV to global vector if it isn't already
  if (base.blacs(ICTXT=a@CTXT)$NPCOL > 1L){
    c <- new("ddmatrix", Data=matrix(out$IPIV, nrow=1L),
              dim=c(1L, b@dim[1L]), ldim=c(1L, b@ldim[1L]), 
              bldim=b@bldim, CTXT=oldctxt)
    pivot <- as.vector(c)
  } else {
    pivot <- out$IPIV
  }
  
  if (out$RANK < n){
#    vec <- as.ddmatrix(matrix(NA, nrow=1, ncol=nrhs), bldim=b@bldim)
    if (m >= n)
      b[(out$RANK+1L):n, , ICTXT=b@CTXT] <- as.double(NA)
    else {
      if (out$RANK < m)
        b[(out$RANK+1L):m, , ICTXT=b@CTXT] <- as.double(NA)
    }
    if (any(pivot - temp != 0L)){
      perm <- sapply(temp, function(i) temp[which(i==pivot)])
      b <- base.redistribute(dx=b, bldim=b@bldim, ICTXT=2L)
      b <- b[perm, , ICTXT=oldctxt]
    } else {
      b <- base.redistribute(dx=b, bldim=b@bldim, ICTXT=oldctxt)
    }
  } else {
    b <- base.redistribute(dx=b, bldim=b@bldim, ICTXT=oldctxt)
  }
  
  # rownames
#  if (base.ownany(dim=b@dim, bldim=b@bldim, CTXT=b@CTXT)){
#    coords <- sapply(temp, function(i) base.g2l_coord(ind=i, dim=b@dim, bldim=b@bldim, ICTXT=b@CTXT)[5])
#    mycoords <- coords[which(!is.na(coords))]
#    
#    rownames(b@Data) <- paste("x", mycoords, sep="")
#  }
  
  qr <- list(qr=a, tau=out$TAU, pivot=pivot, tol=tol, rank=out$RANK)
  
  ret <- list(coefficients=b, residuals=residuals, effects=NULL, 
              rank=out$RANK, fitted.values=fitted.values, assign=NULL,
              qr=qr, df.residual=(a@dim[1] - out$RANK))

  return( ret )
}


# qr.Q()
# recover Q from base.rpdgeqrf
base.pdorgqr <- function(qr)
{
  x <- qr$qr
  
  if (x@dim[1] < x@dim[2])
    x <- x[, 1:qr$rank]
  
  # Matrix descriptors
  desca <- base.descinit(x@dim, x@bldim, x@ldim, ICTXT=x@CTXT)

  m <- desca[3]
  n <- desca[4]
  
  k <- qr$rank
  
  out <- .Call("R_PDORGQR",
            as.integer(m), as.integer(n), as.integer(k),
            x@Data, as.integer(x@ldim), as.integer(desca),
            as.double(qr$tau),
            PACKAGE="pbdBASE")
  
  if (out$INFO!=0)
    warning(paste("ScaLAPACK returned INFO=", out$INFO, "; returned solution is likely invalid", sep=""))

  x@Data <- out$A
  
  return( x )
}


# qr.R
base.qr.R <- function(qr, complete=FALSE)
{
  ret <- qr$qr
  
  if (!complete)
    if (min(ret@dim)!=ret@dim[1])
      ret <- ret[1:min(ret@dim), ]
  
  ret@Data <- base.low2zero(A=ret@Data, dim=ret@dim, ldim=ret@ldim, bldim=ret@bldim, CTXT=ret@CTXT)
  
  # not particularly efficient, but no one should really be calling this...
  rank <- qr$rank
  n <- ret@dim[1]
  p <- ret@dim[2]
  mn <- min(ret@dim)
  if (rank < p){
    if (n>p)
      for (i in (rank+1):mn)
        ret[i,i] <- 0
  }
  
  return(ret)
}


# multiply Q/Q^T against y
base.pdormqr <- function(qr, y, side='L', trans='T')
{
#  x <- base.pdorgqr(qr)
  x <- qr$qr

  # Matrix descriptors
  desca <- base.descinit(x@dim, x@bldim, x@ldim, ICTXT=x@CTXT)
  descb <- base.descinit(y@dim, y@bldim, y@ldim, ICTXT=y@CTXT)

  m <- desca[3]
  n <- y@dim[2]
  k <- qr$rank
  
  IJ <- 1 # IA, JA, IC, JC
  
#  # Determine size of work array 
#  lwork <- .Fortran("PDORMQR",
#            SIDE=as.character(side), TRANS=as.character(trans),
#            M=as.integer(m), N=as.integer(n), K=as.integer(k),
#            A=double(1), as.integer(IJ), as.integer(IJ), DESCA=as.integer(desca), 
#            TAU=double(1), 
#            C=double(1), as.integer(IJ), as.integer(IJ), DESCC=as.integer(descc), 
#            WORK=double(1), LWORK=as.integer(-1), INFO=integer(1),
#            PACKAGE="pbdBASE")$WORK[1]

#  # perform QR
#  out <- .Fortran("PDORMQR",
#            SIDE=as.character(side), TRANS=as.character(trans),
#            M=as.integer(m), N=as.integer(n), K=as.integer(k),
#            A=x@Data, as.integer(IJ), as.integer(IJ), DESCA=as.integer(desca), 
#            TAU=as.double(qr$tau), 
#            C=y@Data, as.integer(IJ), as.integer(IJ), DESCC=as.integer(descc),
#            WORK=double(lwork), LWORK=as.integer(lwork), INFO=integer(1),
#            PACKAGE="pbdBASE")

  out <- .Call("R_PDORMQR",
            as.character(side), as.character(trans),
            as.integer(m), as.integer(n), as.integer(k),
            x@Data, as.integer(x@ldim), as.integer(desca),
            as.double(qr$tau),
            y@Data, as.integer(y@ldim), as.integer(descb),
            PACKAGE="pbdBASE")
  
  
  if (out$INFO!=0)
    warning(paste("ScaLAPACK returned INFO=", out$INFO, "; returned solution is likely invalid", sep=""))

  y@Data <- out$B
  
  return( y )
}



# reduces upper trapezoidal to traingular form
base.pdtzrzf <- function(x)
{
  # Matrix descriptors
  desca <- base.descinit(x@dim, x@bldim, x@ldim, ICTXT=x@CTXT)

  m <- desca[3]
  n <- desca[4]
  
  k <- qr$rank
  
  out <- .Call("R_PDTZRZF",
            as.integer(m), as.integer(n),
            x@Data, as.integer(x@ldim), as.integer(desca),
            as.double(qr$tau),
            PACKAGE="pbdBASE")
  
  if (out$INFO!=0)
    warning(paste("ScaLAPACK returned INFO=", out$INFO, "; returned solution is likely invalid", sep=""))
  
  dim <- rep(x@dim[1L], 2)
  
  ldim <- base.numroc(dim=dim, bldim=x@bldim, ICTXT=x@CTXT, fixme=FALSE)
  if (any(ldim<1))
    Data <- matrix(0)
  else
    Data <- matrix(out$A, nrow=ldim[1], ncol=ldim[2])
  
  ret <- new("ddmatrix", Data=Data, dim=dim,
                    ldim=dim(Data), bldim=x@bldim, CTXT=x@CTXT)
  
  return( ret )
}


# triangle system solve --- probably not needed
base.pdtrsv <- function(x, y, uplo='U', trans='T')
{
  # Matrix descriptors
  desca <- base.descinit(x@dim, x@bldim, x@ldim, ICTXT=x@CTXT)
  descb <- base.descinit(y@dim, y@bldim, y@ldim, ICTXT=y@CTXT)
  
  n <- desca[4]
  
  dg <- 'N'
  
#pdtrsv(uplo, trans, diag, n, a, ia, ja, desca, x, ix, jx, descx, incx)

  # perform QR
#  out <- .C("pdtrsv_",
#            UPLO=as.character(uplo), TRANS=as.character(trans),
#            DIAG=as.character('N'), N=as.integer(n),
#            A=x@Data, as.integer(1), as.integer(1), DESCA=as.integer(desca), 
#            X=y@Data, as.integer(1), as.integer(1), DESCC=as.integer(descc),
#            INCX=as.integer(1),
#            PACKAGE="pbdBASE")

  out <- .Call("R_PDTRSV",
            as.character(uplo), as.character(trans), as.character(dg),
            as.integer(n),
            x@Data, as.integer(x@ldim), as.integer(desca),
            y@Data, as.integer(y@ldim), as.integer(descb),
            PACKAGE="pbdBASE")
  
  if (out$INFO!=0)
    warning(paste("ScaLAPACK returned INFO=", out$INFO, "; returned solution is likely invalid", sep=""))

  y@Data <- out$C
  
  return( y )
}

