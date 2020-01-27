# ------------------------------------------------
# Linear Equations
# ------------------------------------------------

#' rpdgetri
#' 
#' Matrix inversion.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param n
#' Problem size.
#' @param a
#'  Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' 
#' @useDynLib pbdBASE R_PDGETRI
#' @export
base.rpdgetri <- function(n, a, desca)
{
  aldim <- base.numroc(desca[3:4], desca[5:6], ICTXT=desca[2])
  
  if (!is.double(a))
    storage.mode(a) <- "double"
  
  out <- .Call(R_PDGETRI, a, as.integer(desca))
  
  if (out$info!=0)
    pbdMPI::comm.warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))
  
  out$A
}



#' rpdgesv
#' 
#' Solving a (square) system of equations.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param n
#' Problem size.
#' @param nrhs
#' Number of right hand sides.
#' @param a,b
#' Matrix.
#' @param desca,descb
#' ScaLAPACK descriptor array.
#' 
#' @useDynLib pbdBASE R_PDGESV
#' @export
base.rpdgesv <- function(n, nrhs, a, desca, b, descb)
{
  aldim <- base.numroc(desca[3:4], desca[5:6], ICTXT=desca[2])
  bldim <- base.numroc(descb[3:4], descb[5:6], ICTXT=descb[2])
  
  # max of the local dimensions
  mxldims <- c(base.maxdim(aldim), base.maxdim(bldim))
  
  if (!is.double(a))
    storage.mode(a) <- "double"
  if (!is.double(b))
    storage.mode(b) <- "double"
  
  # Call ScaLAPACK
  out <- .Call(R_PDGESV,
               as.integer(n), as.integer(nrhs), as.integer(mxldims),
               a, as.integer(desca), b, as.integer(descb))
  
  if (out$info!=0)
    pbdMPI::comm.warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))
  
  out$B 
}



# ------------------------------------------------
# Matrix Factorizations
# ------------------------------------------------

#' rpdgesvd
#' 
#' SVD.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param jobu,jobvt
#' Control for u/vt return.
#' @param m,n
#' Problem size.
#' @param a
#'  Matrix.
#' @param desca,descu,descvt
#' ScaLAPACK descriptor array.
#' @param ...
#' Ignored
#' @param inplace
#' Should the computation be done in-place or not.  For REALLY advanced users only.
#' @param comm
#' An MPI (not BLACS) communicator.
#' 
#' @useDynLib pbdBASE R_PDGESVD
#' @export
base.rpdgesvd <- function(jobu, jobvt, m, n, a, desca, descu, descvt, ..., inplace=FALSE, comm = .pbd_env$SPMD.CT$comm)
{
  size <- min(m, n)
  
  aldim <- dim(a)
  uldim <- base.numroc(descu[3:4], descu[5:6], ICTXT=descu[2])
  vtldim <- base.numroc(descvt[3:4], descvt[5:6], ICTXT=descvt[2])
  
  mxa <- pbdMPI::allreduce(max(aldim), op='max', comm=comm)
  mxu <- pbdMPI::allreduce(max(uldim), op='max', comm=comm)
  mxvt <- pbdMPI::allreduce(max(vtldim), op='max', comm=comm)
  
  if (all(aldim==1))
    desca[9L] <- mxa
  if (all(uldim==1))
    descu[9L] <- mxu
  if (all(vtldim==1))
    descvt[9L] <- mxvt
  
  if (desca[3L]>1){
    if (pbdMPI::allreduce(desca[9L], op='max', comm=comm)==1)
      desca[9L] <- mxa
  }
  if (descu[3L]>1){
    if (pbdMPI::allreduce(descu[9L], op='max', comm=comm)==1)
      desca[9L] <- mxu
  }
  if (descvt[3L]>1){
    if (pbdMPI::allreduce(descvt[9L], op='max', comm=comm)==1)
      desca[9L] <- mxvt
  }
  
  if (!is.double(a))
    storage.mode(a) <- "double"
  
  # FIXME currently does nothing
  if (inplace)
    inplace <- 'Y'
  else
    inplace <- 'N'
  
  # Call ScaLAPACK
  out <- .Call(R_PDGESVD, 
              as.integer(m), as.integer(n), as.integer(size),
              a, as.integer(desca), 
              as.integer(uldim), as.integer(descu),
              as.integer(vtldim), as.integer(descvt),
              as.character(jobu), as.character(jobvt), inplace)
  
  if (out$info!=0)
    pbdMPI::comm.warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""), comm=comm)
  
  ret <- list( d=out$d, u=out$u, vt=out$vt )
  ret
}



#' rpdsyevr
#' 
#' Symmetric eigenvalue decomposition.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param jobz
#' Control for if vectors/values/both are returned.
#' @param uplo
#' Triangle where the information is stored (in the symmetric matrix).
#' @param n
#' Problem size.
#' @param a
#' Matrix.
#' @param desca,descz
#' ScaLAPACK descriptor array.
#' 
#' @useDynLib pbdBASE R_PDSYEVR
#' @export
base.rpdsyevr <- function(jobz, uplo, n, a, desca, descz)
{
  aldim <- dim(a)
  # zldim <- base.numroc(descz[3:4], descz[5:6], ICTXT=descz[2])
  
  # mxa <- pbdMPI::allreduce(max(aldim), op='max')
  # mxz <- pbdMPI::allreduce(max(zldim), op='max')
  # 
  # if (all(aldim==1))
  #     desca[9L] <- mxa
  # if (all(zldim==1))
  #     descz[9L] <- mxz
  
  if (!is.double(a))
    storage.mode(a) <- "double"
  
  # Call ScaLAPACK
  out <- .Call(R_PDSYEVR, as.character(jobz), as.character(uplo),
      as.integer(n), a, as.integer(desca), as.integer(descz))
  
  if (out$info!=0)
    pbdMPI::comm.warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))
  
  out$info <- NULL
  out
}



#' rpdpotrf
#' 
#' Cholesky factorization.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param uplo
#' Triangle where the information is stored (in the symmetric matrix).
#' @param n
#' Problem size.
#' @param a
#'  Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' 
#' @useDynLib pbdBASE R_PDPOTRF
#' @export
base.rpdpotrf <- function(uplo, n, a, desca)
{
  if (!is.double(a))
    storage.mode(a) <- "double"
  
  # Call ScaLAPACK
  ret <- .Call(R_PDPOTRF,
               as.integer(n), a, as.integer(desca),
               as.character(uplo))
  
  if (ret$info!=0)
    pbdMPI::comm.warning(paste("ScaLAPACK returned INFO=", ret$info, "; returned solution is likely invalid", sep=""))
  
  ret
}



#' rpdsyevx
#' 
#' Genearlized eigenvalue problem.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param jobz
#' Control for if vectors/values/both are returned.
#' @param range
#' Parameter to determine the search criteria for eigenvalues.
#' @param n
#' Problem size.
#' @param a
#'  Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' @param vl,vu
#' Endpoints of the interval subset of the real line in which to search for eigenvalues, if specified by \code{range}.
#' @param il,iu
#' Eigenvalues with indices \code{il}, ..., \code{iu} will be found, if specified by \code{range}.
#' @param abstol
#' Absolute error tolerance for the eigenvalues.
#' @param orfac
#' Eigenvectors with eigenvalues below orfac*norm(a) of each other are reorthogonalized.
#' 
#' @useDynLib pbdBASE R_PDSYEVX
#' @export
base.rpdsyevx <- function(jobz, range, n, a, desca, vl, vu, il, iu, abstol=1e-8, orfac=1e-3)
{
  if (!is.double(a))
    storage.mode(a) <- "double"
  
  
  # Call ScaLAPACK
  ret <- .Call(R_PDSYEVX,
             as.character(jobz), as.character(range),
             as.integer(n), a, as.integer(desca),
             as.double(vl), as.double(vu), as.integer(il), as.integer(iu),
             as.double(abstol), as.double(orfac))
  
  ret
}


#' rpdgetrf
#' 
#' LU factorization.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param a
#'  Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' 
#' @useDynLib pbdBASE R_PDGETRF
#' @export
base.rpdgetrf <- function(a, desca)
{
  m <- desca[3L]
  n <- desca[4L]
  
  aldim <- dim(a)
  
  mxrow <- base.igamx2d(ICTXT=desca[2L], SCOPE='Row', m=1L, n=1L, x=aldim[1L], lda=1L, RDEST=-1L, CDEST=-1L)
  lipiv <- mxrow + desca[5L]
#  lipiv <- base.maxdim(aldim)[1L] + desca[5L]
  
  if (!is.double(a))
    storage.mode(a) <- "double"
  
  # Call ScaLAPACK
  out <- .Call(R_PDGETRF,
               as.integer(m), as.integer(n),
               a, as.integer(aldim), as.integer(desca),
               as.integer(lipiv))
  
  if (out$info!=0)
    pbdMPI::comm.warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))
  
  out$A
}



# ------------------------------------------------
# Auxillary
# ------------------------------------------------

#' rpdlange
#' 
#' Matrix norms.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param norm
#' Type of norm.
#' @param m,n
#' Problem size
#' @param a
#' Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' 
#' @useDynLib pbdBASE R_PDLANGE
#' @export
base.rpdlange <- function(norm, m, n, a, desca)
{
  if (length(norm)>1L)
    norm <- norm[1L]
  
  norm <- toupper(norm)
  
  if (!is.double(a))
    storage.mode(a) <- "double"
  
  ret <- .Call(R_PDLANGE, 
              norm, as.integer(m), as.integer(n),
              a, as.integer(desca))
  
  ret
}



#' rpdtrcon
#' 
#' Inverse condition number of a triangular matrix.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param norm
#' Type of norm.
#' @param uplo
#' Triangle where information is stored.
#' @param diag
#' Specifies if the matrix is unit triangular or not.
#' @param n
#' Problem size
#' @param a
#' Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' 
#' @useDynLib pbdBASE R_PDTRCON
#' @export
base.rpdtrcon <- function(norm, uplo, diag, n, a, desca)
{
  if (length(norm)>1L)
    norm <- norm[1L]
  
  norm <- toupper(norm)
  uplo <- toupper(uplo)
  diag <- toupper(diag)
  
  if (!is.double(a))
    storage.mode(a) <- "double"
  
  ret <- .Call(R_PDTRCON, 
              norm, uplo, diag, 
              as.integer(n), a, as.integer(desca))
  
  if (ret[2L] < 0)
    pbdMPI::comm.warning(paste("INFO =", ret[2L]))
  
  ret[1L]
}



#' rpdgecon
#' 
#' Inverse condition number of a general matrix.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param norm
#' Type of norm.
#' @param m,n
#' Problem size
#' @param a
#' Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' 
#' @useDynLib pbdBASE R_PDGECON
#' @export
base.rpdgecon <- function(norm, m, n, a, desca)
{
  if (length(norm)>1L)
    norm <- norm[1L]
  
  norm <- toupper(norm)
  
  if (!is.double(a))
    storage.mode(a) <- "double"
  
  ret <- .Call(R_PDGECON, norm, as.integer(m), as.integer(n), a, as.integer(desca))
  
  if (ret[2] < 0)
    pbdMPI::comm.warning(paste("INFO =", ret[2]))
  
  ret[1]
}



#' det
#' 
#' Determinant.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param a
#' Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' 
#' @useDynLib pbdBASE R_det
#' @export
base.det = function(a, desca)
{
  if (!is.double(a))
    storage.mode(a) = "double"
  
  ret = .Call(R_det, a, as.integer(desca))
  
  if (ret$info!=0)
    pbdMPI::comm.warning(paste("ScaLAPACK returned INFO=", ret$info, "; returned solution is likely invalid", sep=""))
  
  ret$info = NULL
  ret
}
