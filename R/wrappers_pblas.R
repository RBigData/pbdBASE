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

