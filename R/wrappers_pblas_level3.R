# ################################################
# ------------------------------------------------
# PDTRAN:  Matrix transpose
# ------------------------------------------------
# ################################################

base.pdtran <- function(a)
{
  ICTXT <- a@CTXT
  blacs_ <- base.blacs(ICTXT=ICTXT)
  MYROW <- blacs_$MYROW
  MYCOL <- blacs_$MYCOL
  
  m <- a@dim[2]
  n <- a@dim[1]
  
  desca <- base.descinit(a@dim, a@bldim, a@ldim, ICTXT=ICTXT)

  cdim <- c(m, n)
  cldim <- base.numroc(cdim, a@bldim, ICTXT=ICTXT)

  c <- new("ddmatrix", Data=matrix(nrow=0, ncol=0),
                       dim=cdim, ldim=cldim, bldim=a@bldim, CTXT=ICTXT)

  descc <- base.descinit(c@dim, c@bldim, c@ldim, ICTXT=ICTXT)

  ret <- .Call("R_PDTRAN",
                as.integer(m), as.integer(n),
                a@Data, as.integer(desca),
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

base.pdgemm <- function(a, b)
{
  ICTXT <- a@CTXT
  
  m <- a@dim[1L]
  n <- b@dim[2L]
  k <- b@dim[1L]
  
  bldim <- a@bldim
  
  cdim <- c(a@dim[1L], b@dim[2L])
  cldim <- base.numroc(cdim, a@bldim, ICTXT=ICTXT)
  
  desca <- base.descinit(a@dim, bldim, a@ldim, ICTXT=ICTXT)
  descb <- base.descinit(b@dim, bldim, b@ldim, ICTXT=ICTXT)
  descc <- base.descinit(cdim, bldim, cldim, ICTXT=ICTXT)
  
  trans <- 'N'
  
  ret <- .Call("R_PDGEMM",
                  as.character(trans), as.character(trans),
                  as.integer(m), as.integer(n), as.integer(k),
                  a@Data, as.integer(desca),
                  b@Data, as.integer(descb),
                  as.integer(cldim), as.integer(descc),
                  PACKAGE="pbdBASE"
                 )
  
  c <- new("ddmatrix", Data=ret, dim=cdim, ldim=cldim, bldim=a@bldim, CTXT=ICTXT)
  
  return(c)
}
