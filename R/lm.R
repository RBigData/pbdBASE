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
  
  # Sometimes R mistakenly frees these matrices...
  if (!base.ownany(dim=desca[5L:6L], bldim=desca[5L:6L], ICTXT=desca[2L]))
    ret$A <- matrix(0.0, 1L, 1L)
  
  if (!base.ownany(dim=descb[5L:6L], bldim=descb[5L:6L], ICTXT=descb[2L]))
    ret$EFF <- ret$RSD <- ret$FT <- matrix(0.0, 1L, 1L)
  
  
  if (ret$INFO!=0)
    comm.warning(paste("ScaLAPACK returned INFO=", ret$INFO, "; returned solution is likely invalid", sep=""))
  
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
  
  if (comm.rank()!=0)
    rank <- 0L
  else
    rank <- ret$rank
    
  rank <- pbdMPI::allreduce(rank)
  
  if (ret$INFO!=0)
    comm.warning(paste("ScaLAPACK returned INFO=", ret$INFO, "; returned solution is likely invalid", sep=""))
  
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
    comm.warning(paste("ScaLAPACK returned INFO=", out$INFO, "; returned solution is likely invalid", sep=""))
  
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
    comm.warning(paste("ScaLAPACK returned INFO=", out$INFO, "; returned solution is likely invalid", sep=""))
  
  return( out$B )
}

