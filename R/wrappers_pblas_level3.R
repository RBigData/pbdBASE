# ################################################
# ------------------------------------------------
# Level 3 PBLAS
# ------------------------------------------------
# ################################################


# ------------------------------------------------
# PDTRAN:  Matrix transpose
# ------------------------------------------------

base.rpdtran <- function(x)
{
  ICTXT <- x@CTXT
  blacs_ <- base.blacs(ICTXT=ICTXT)
  
  m <- x@dim[2]
  n <- x@dim[1]
  
  desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=ICTXT)

  cdim <- c(m, n)
  cldim <- base.numroc(cdim, x@bldim, ICTXT=ICTXT)

  c <- new("ddmatrix", Data=matrix(nrow=0, ncol=0),
                       dim=cdim, ldim=cldim, bldim=x@bldim, CTXT=ICTXT)

  descc <- base.descinit(c@dim, c@bldim, c@ldim, ICTXT=ICTXT)

  if (!is.double(x@Data))
    storage.mode(x@Data) <- "double"

  ret <- .Call("R_PDTRAN",
                as.integer(m), as.integer(n),
                x@Data, as.integer(desca),
                as.integer(cldim), as.integer(descc),
                PACKAGE="pbdBASE"
              )

  c@Data <- ret

  return(c)
}

# ------------------------------------------------
# PDGEMM:  Matrix-Matrix multiplication
# ------------------------------------------------

base.rpdgemm <- function(transx='N', transy='N', x, y, outbldim=x@bldim)
{
  if (length(outbldim)==1L)
    outbldim <- rep(outbldim, 2)
  
  ICTXT <- x@CTXT
  
  m <- x@dim[1L]
  n <- y@dim[2L]
  k <- y@dim[1L]
  
  bldimx <- x@bldim
  bldimy <- y@bldim
  
  cdim <- c(x@dim[1L], y@dim[2L])
  cldim <- base.numroc(cdim, outbldim, ICTXT=ICTXT)
  
  desca <- base.descinit(dim=x@dim, bldim=bldimx, ldim=x@ldim, ICTXT=ICTXT)
  descb <- base.descinit(dim=y@dim, bldim=bldimy, ldim=y@ldim, ICTXT=ICTXT)
  descc <- base.descinit(dim=cdim, bldim=outbldim, ldim=cldim, ICTXT=ICTXT)
  
  if (!is.double(x@Data))
    storage.mode(x@Data) <- "double"
  if (!is.double(y@Data))
    storage.mode(y@Data) <- "double"
    
  ret <- .Call("R_PDGEMM",
                  as.character(transx), as.character(transy),
                  as.integer(m), as.integer(n), as.integer(k),
                  x@Data, as.integer(desca),
                  y@Data, as.integer(descb),
                  as.integer(cldim), as.integer(descc),
                  PACKAGE="pbdBASE"
                 )
  
  c <- new("ddmatrix", Data=ret, dim=cdim, ldim=cldim, bldim=outbldim, CTXT=ICTXT)
  
  return(c)
}

# ------------------------------------------------
# PDSVRK:  Symmetric Rank-k Update
# ------------------------------------------------

base.rpdsvrk <- function(trans, x, outbldim=x@bldim)
{
  if (length(outbldim)==1L)
    outbldim <- rep(outbldim, 2)
  
  ICTXT <- x@CTXT
  
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
  
  c <- new("ddmatrix", Data=ret, dim=cdim, ldim=cldim, bldim=outbldim, CTXT=ICTXT)
  
  return(c)
}


base.crossprod <- function(trans, x)
{
  ICTXT <- x@CTXT
  
  if (trans=='N' || trans=='n'){
    n <- x@dim[2L]
    k <- x@dim[1L]
  } else {
    n <- x@dim[1L]
    k <- x@dim[2L]
  }
  
  bldim <- x@bldim
  
  cdim <- c(n, n)
  cldim <- base.numroc(cdim, bldim, ICTXT=ICTXT)
  
  desca <- base.descinit(dim=x@dim, bldim=bldim, ldim=x@ldim, ICTXT=ICTXT)
  descc <- base.descinit(dim=cdim, bldim=bldim, ldim=cldim, ICTXT=ICTXT)
  
  if (!is.double(x@Data))
    storage.mode(x@Data) <- "double"
  
  ret <- .Call("R_PDCROSSPROD",
                  as.character(trans),
                  x@Data, as.integer(desca),
                  as.integer(cldim), as.integer(descc),
                  PACKAGE="pbdBASE"
                 )
  
  c <- new("ddmatrix", Data=ret, dim=cdim, ldim=cldim, bldim=bldim, CTXT=ICTXT)
  
  return(c)
}

