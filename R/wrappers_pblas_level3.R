# ################################################
# ------------------------------------------------
# PDTRAN:  Matrix transpose
# ------------------------------------------------
# ################################################

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

  ret <- .Call("R_PDTRAN",
                as.integer(m), as.integer(n),
                x@Data, as.integer(desca),
                as.integer(cldim), as.integer(descc),
                PACKAGE="pbdBASE"
              )

  c@Data <- ret

  return(c)
}

# ################################################
# ------------------------------------------------
# PDGEMM:  Matrix-Matrix multiplication
# ------------------------------------------------
# ################################################

base.rpdgemm <- function(x, y, outbldim=x@bldim)
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
  
  trans <- 'N'
  
  ret <- .Call("R_PDGEMM",
                  as.character(trans), as.character(trans),
                  as.integer(m), as.integer(n), as.integer(k),
                  x@Data, as.integer(desca),
                  y@Data, as.integer(descb),
                  as.integer(cldim), as.integer(descc),
                  PACKAGE="pbdBASE"
                 )
  
  c <- new("ddmatrix", Data=ret, dim=cdim, ldim=cldim, bldim=outbldim, CTXT=ICTXT)
  
  return(c)
}
