# Work in progress --- not ready for the big time

# Wrapper for custom PDGELS function, which solves linear least
# squares problems.
base.rpdgels <- function(a, b, tol=1e-7)
{
  # BLACS stuff
  ICTXT <- a@CTXT

  # Matrix descriptors
  desca <- base.descinit(a@dim, a@bldim, a@ldim, ICTXT=ICTXT)
  descb <- base.descinit(b@dim, b@bldim, b@ldim, ICTXT=ICTXT)

  m <- desca[3]
  n <- desca[4]
  nrhs <- descb[4]
  
  # Leave this alone

  # Determine size of work array 
  lwork <- .Fortran("RPDGELS", 
                    as.double(tol), as.character("N"), 
                    as.integer(m), as.integer(n), as.integer(nrhs),
                    double(1), as.integer(1), as.integer(1), as.integer(desca),
                    double(1), as.integer(1), as.integer(1), as.integer(descb),
                    WORK=double(1), as.integer(-1), as.integer(1), 
                    integer(1), INFO=as.integer(0)
                    )$WORK[1]

  # Convert to .Call()

  # RPDGELS( TOL, TRANS, M, N, NRHS, A, IA, JA, DESCA, 
  #          B, IB, JB, DESCB, WORK, LWORK, IPIV, RANK, INFO )
  
  # in: tol, trans, m, n, nrhs, IA, JA, desca, IB, JB, descb, lwork
  # in/out: A, B, 
  # out: IPIV, RANK
  # local (just allocate in C, not a SEXP): work

  # Fit the model
  out <- .Fortran("RPDGELS", 
                  as.double(tol), as.character("N"),
                  as.integer(m), as.integer(n), as.integer(nrhs),
                  A=a@Data, as.integer(1), as.integer(1), as.integer(desca), 
                  B=b@Data, as.integer(1), as.integer(1), as.integer(descb),
                  double(lwork), as.integer(lwork), IPIV=integer(a@ldim[2]),
                  RANK=integer(1), INFO=as.integer(0),
                  PACKAGE="pbdBASE")

  if (out$INFO!=0)
    warning(paste("ScaLAPACK returned INFO=", out$INFO, "; returned solution is likely invalid", sep=""))

  return( list(QR=out$A, X=out$B, IPIV=out$IPIV) )
}


