# ################################################
# ------------------------------------------------
# Level 3 PBLAS
# ------------------------------------------------
# ################################################


# ------------------------------------------------
# PDTRAN:  Matrix transpose
# ------------------------------------------------

base.rpdtran <- function(m, n, a, desca, descc)
{
  if (!is.double(a))
    storage.mode(a) <- "double"
  
  cldim <- base.numroc(descc[3:4], descc[5:6], ICTXT=descc[2])
  
  ret <- .Call("R_PDTRAN",
                as.integer(m), as.integer(n),
                a, as.integer(desca),
                as.integer(cldim), as.integer(descc),
                PACKAGE="pbdBASE"
              )
  
  return(ret)
}

# ------------------------------------------------
# PDGEMM:  Matrix-Matrix multiplication
# ------------------------------------------------

base.rpdgemm <- function(transx, transy, m, n, k, x, descx, y, descy, descc)
{
  transx <- toupper(transx)
  transy <- toupper(transy)
  
  cldim <- base.numroc(descc[3:4], descc[5:6], ICTXT=descc[2])
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
    
  ret <- .Call("R_PDGEMM",
                transx, transy,
                as.integer(m), as.integer(n), as.integer(k),
                x, as.integer(descx),
                y, as.integer(descy),
                as.integer(cldim), as.integer(descc),
                PACKAGE="pbdBASE")
  
  return( ret )
}

# ------------------------------------------------
# PDSVRK:  Symmetric Rank-k Update
# ------------------------------------------------

base.rpdsvrk <- function(trans, x, outbldim=x@bldim)
{
  if (length(outbldim)==1L)
    outbldim <- rep(outbldim, 2)
  
  ICTXT <- x@ICTXT
  
  if (trans=='N' || trans=='n'){
    n <- x@dim[1L]
    k <- x@dim[2L]
  } else {
    n <- x@dim[2L]
    k <- x@dim[1L]
  }
  
  bldim <- x@bldim
  
  cdim <- c(n, n)
  cldim <- base.numroc(cdim, outbldim, ICTXT=ICTXT)
  
  desca <- base.descinit(dim=x@dim, bldim=bldim, ldim=x@ldim, ICTXT=ICTXT)
  descc <- base.descinit(dim=cdim, bldim=outbldim, ldim=cldim, ICTXT=ICTXT)
  
  if (!is.double(x@Data))
    storage.mode(x@Data) <- "double"
  
  uplo <- 'U'
  
  ret <- .Call("R_PDSYRK",
                  as.character(uplo), as.character(trans),
                  as.integer(n), as.integer(k),
                  x@Data, as.integer(desca),
                  as.integer(cldim), as.integer(descc),
                  PACKAGE="pbdBASE"
                 )
  
  c <- new("ddmatrix", Data=ret, dim=cdim, ldim=cldim, bldim=outbldim, ICTXT=ICTXT)
  
  return(c)
}


base.crossprod <- function(trans, x, descx, descc)
{
  trans <- toupper(trans)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  cldim <- base.numroc(descc[3:4], descc[5:6], ICTXT=descc[2])
  
  ret <- .Call("R_PDCROSSPROD",
                  trans, x, as.integer(descx),
                  as.integer(cldim), as.integer(descc),
                  PACKAGE="pbdBASE")
  
  return( ret )
}

