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
                  a@Data, as.integer(cldim),
                  as.integer(ICTXT), as.integer(MYROW), as.integer(MYCOL),
                  as.integer(desca), as.integer(descc),
                  as.integer(m), as.integer(n),
                  PACKAGE="pbdBASE"
                 )

#  ret <- .Call("R_PDTRAN",
#                  as.integer(m), as.integer(n),
#                  a@Data, as.integer(desca), 
#                  as.integer(cldim), as.integer(descc),
#                  as.integer(ICTXT),
#                  PACKAGE="pbdBASE"
#                 )

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
  blacs_ <- base.blacs(ICTXT=a@CTXT)
  MYROW <- blacs_$MYROW
  MYCOL <- blacs_$MYCOL

  m <- a@dim[1]
  n <- b@dim[2]
  k <- b@dim[1]

  cdim <- c(a@dim[1], b@dim[2])
  
  cldim <- base.numroc(cdim, a@bldim, ICTXT=ICTXT)
  
  c <- new("ddmatrix", Data=matrix(nrow=0, ncol=0),
                       dim=cdim, ldim=cldim, bldim=a@bldim, CTXT=ICTXT)

  desca <- base.descinit(a@dim, a@bldim, a@ldim, ICTXT=ICTXT)
  descb <- base.descinit(b@dim, b@bldim, b@ldim, ICTXT=ICTXT)
  descc <- base.descinit(c@dim, c@bldim, c@ldim, ICTXT=ICTXT)
  
  c@Data <- .Call("R_PDGEMM",
                  a@Data, b@Data, as.integer(cldim),
                  as.integer(ICTXT), as.integer(MYROW), as.integer(MYCOL),
                  as.integer(desca), as.integer(descb), as.integer(descc),
                  as.integer(m), as.integer(n), as.integer(k),
                  PACKAGE="pbdBASE"
                 )
  return(c)
}
