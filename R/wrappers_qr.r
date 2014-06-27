# qr()
base.rpdgeqpf <- function(tol, m, n, x, descx)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  ret <- .Call("R_PDGEQPF",
               as.double(tol), as.integer(m), as.integer(n), x, as.integer(descx),
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

