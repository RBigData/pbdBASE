# lm.fit()
# Wrapper for custom PDGELS function, which solves linear least
# squares problems.
base.rpdgels <- function(tol, m, n, nrhs, a, desca, b, descb)
{
  # FIXME adjustment for weird lda issue
  mxlda <- pbdMPI::allreduce(desca[9], op='max')
  mxldb <- pbdMPI::allreduce(descb[9], op='max')
  
  if (desca[9]==1)
    desca[9] <- mxlda
  if (descb[9]==1)
    descb[9] <- mxldb
  
  if (!is.double(a))
    storage.mode(a) <- "double"
  if (!is.double(b))
    storage.mode(b) <- "double"
  
  aldim <- as.integer(dim(a))
  bldim <- as.integer(dim(b))
  
  ltau <- as.integer(min(m, n))
  
  ret <- .Call("R_PDGELS",
            TOL=as.double(tol), M=as.integer(m), N=as.integer(n), NRHS=as.integer(nrhs),
            A=a, ALDIM=as.integer(aldim), DESCA=as.integer(desca),
            B=b, BLDIM=as.integer(bldim), DESCB=as.integer(descb),
            LTAU=ltau,
            PACKAGE="pbdBASE")
  
  if (ret$INFO!=0)
    warning(paste("ScaLAPACK returned INFO=", out$INFO, "; returned solution is likely invalid", sep=""))
  
  
  return( ret )
}



# qr()
base.rpdgeqpf <- function(tol, m, n, x, descx)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  ret <- .Call("R_PDGEQPF",
               as.double(tol), as.integer(m), as.integer(n), 
               x, as.integer(dim(x)), as.integer(descx),
               as.integer(min(m, n)),
               PACKAGE="pbdBASE")
  
  if (ret$INFO!=0)
    warning(paste("ScaLAPACK returned INFO=", ret$INFO, "; returned solution is likely invalid", sep=""))
  
  
  return( ret )
}




# qr.Q()
# recover Q from base.rpdgeqrf
base.rpdorgqr <- function(m, n, k, qr, descqr, tau)
{
  if (!is.double(qr))
    storage.mode(qr) <- "double"
  
  if (!is.double(tau))
    storage.mode(tau) <- "double"
  
  out <- .Call("R_PDORGQR",
            as.integer(m), as.integer(n), as.integer(k),
            qr, as.integer(dim(qr)), as.integer(descqr), tau,
            PACKAGE="pbdBASE")
  
  if (out$INFO!=0)
    warning(paste("ScaLAPACK returned INFO=", out$INFO, "; returned solution is likely invalid", sep=""))
  
  ret <- out$A
  
  return( ret )
}



# multiply Q/Q^T against y
base.rpdormqr <- function(side, trans, m, n, k, qr, descqr, tau, c, descc)
{
  # FIXME adjustment for weird lda issue
  mxlda <- pbdMPI::allreduce(descqr[9], op='max')
  mxldb <- pbdMPI::allreduce(descc[9], op='max')
  
  if (descqr[9]==1)
    descqr[9] <- mxlda
  if (descc[9]==1)
    descc[9] <- mxldb
  
  if (!is.double(qr))
    storage.mode(qr) <- "double"
  if (!is.double(c))
    storage.mode(c) <- "double"
  if (!is.double(tau))
    storage.mode(tau) <- "double"
  
  out <- .Call("R_PDORMQR",
            as.character(side), as.character(trans),
            as.integer(m), as.integer(n), as.integer(k),
            qr, as.integer(dim(qr)), as.integer(descqr),
            tau,
            c, as.integer(dim(c)), as.integer(descc),
            PACKAGE="pbdBASE")
  
  if (out$INFO!=0)
    warning(paste("ScaLAPACK returned INFO=", out$INFO, "; returned solution is likely invalid", sep=""))
  
  return( out$B )
}



# reduces upper trapezoidal to traingular form
base.rpdtzrzf <- function(qr)
{
  x <- qr$qr
  
  # Matrix descriptors
  desca <- base.descinit(x@dim, x@bldim, x@ldim, ICTXT=x@ICTXT)

  m <- desca[3]
  n <- desca[4]
  
  if (!is.double(x@Data))
    storage.mode(x@Data) <- "double"
  
  out <- .Call("R_PDTZRZF",
            as.integer(m), as.integer(n),
            x@Data, as.integer(x@ldim), as.integer(desca),
            as.double(qr$tau),
            PACKAGE="pbdBASE")
  
  if (out$INFO!=0)
    warning(paste("ScaLAPACK returned INFO=", out$INFO, "; returned solution is likely invalid", sep=""))
  
  dim <- rep(x@dim[1L], 2)
  
  ldim <- base.numroc(dim=dim, bldim=x@bldim, ICTXT=x@ICTXT, fixme=FALSE)
  if (any(ldim<1))
    Data <- matrix(0)
  else
    Data <- matrix(out$A, nrow=ldim[1], ncol=ldim[2])
  
  ret <- new("ddmatrix", Data=Data, dim=dim,
                    ldim=dim(Data), bldim=x@bldim, CTXT=x@ICTXT)
  
  return( ret )
}


# triangle system solve --- probably not needed
base.rpdtrsv <- function(x, y, uplo='U', trans='T')
{
  # Matrix descriptors
  desca <- base.descinit(x@dim, x@bldim, x@ldim, ICTXT=x@ICTXT)
  descb <- base.descinit(y@dim, y@bldim, y@ldim, ICTXT=y@ICTXT)
  
  n <- desca[4]
  
  dg <- 'N'
  
  if (!is.double(x@Data))
    storage.mode(x@Data) <- "double"
  if (!is.double(y@Data))
    storage.mode(y@Data) <- "double"
  
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

