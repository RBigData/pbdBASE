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
  
  eff <- new("ddmatrix", Data=out$EFF, dim=b@dim, 
             ldim=b@bldim, bldim=b@bldim, CTXT=b@CTXT)
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
  
  attr(qr, "class") <- "qr"
  
  ret <- list(coefficients=b, residuals=residuals, effects=eff, 
              rank=out$RANK, fitted.values=fitted.values, assign=attr(a@Data, "assign"),
              qr=qr, df.residual=(a@dim[1] - out$RANK))
  
  return( ret )
}



# qr()
base.rpdgeqpf <- function(x, tol=1e-7)
{
  # Matrix descriptors
  desca <- base.descinit(x@dim, x@bldim, x@ldim, ICTXT=x@CTXT)

  m <- desca[3]
  n <- desca[4]
  
  ret <- .Call("R_PDGEQPF",
               as.double(tol), as.integer(m), as.integer(n), 
               x@Data, as.integer(x@ldim), as.integer(desca),
               as.integer(min(m, n)),
               PACKAGE="pbdBASE")
  
  if (ret$INFO!=0)
    warning(paste("ScaLAPACK returned INFO=", ret$INFO, "; returned solution is likely invalid", sep=""))
  
  x@Data <- ret$qr
  ret$qr <- x
  ret$INFO <- NULL
  
  attr(ret, "class") <- "qr"
  
  return( ret )
}




# qr.Q()
# recover Q from base.rpdgeqrf
base.rpdorgqr <- function(qr)
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
  
  if (!complete){
    if (min(ret@dim)!=ret@dim[1])
      ret <- ret[1:min(ret@dim), ]
  }
  
  ret@Data <- base.low2zero(A=ret@Data, dim=ret@dim, ldim=ret@ldim, bldim=ret@bldim, CTXT=ret@CTXT)
  
  # not particularly efficient, but I don't expect this to get any real use...
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
base.rpdormqr <- function(qr, y, side='L', trans='T')
{
  x <- qr$qr

  # Matrix descriptors
  desca <- base.descinit(x@dim, x@bldim, x@ldim, ICTXT=x@CTXT)
  descb <- base.descinit(y@dim, y@bldim, y@ldim, ICTXT=y@CTXT)

  m <- desca[3]
  n <- y@dim[2]
  k <- qr$rank
  
  # FIXME adjustment for weird lda issue
  mxlda <- pbdMPI::allreduce(desca[9], op='max')
  mxldb <- pbdMPI::allreduce(descb[9], op='max')
  
  if (desca[9]==1)
    desca[9] <- mxlda
  if (descb[9]==1)
    descb[9] <- mxldb
  
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
base.rpdtzrzf <- function(qr)
{
  x <- qr$qr
  
  # Matrix descriptors
  desca <- base.descinit(x@dim, x@bldim, x@ldim, ICTXT=x@CTXT)

  m <- desca[3]
  n <- desca[4]
  
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
base.rpdtrsv <- function(x, y, uplo='U', trans='T')
{
  # Matrix descriptors
  desca <- base.descinit(x@dim, x@bldim, x@ldim, ICTXT=x@CTXT)
  descb <- base.descinit(y@dim, y@bldim, y@ldim, ICTXT=y@CTXT)
  
  n <- desca[4]
  
  dg <- 'N'
  
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

