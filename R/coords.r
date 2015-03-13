### global-to-local
#' @export
indxg2l <- function(INDXGLOB, NB, IPROC, ISRCPROC, NPROCS)
{
  indx <- NB*floor((INDXGLOB - 1L)/(NB*NPROCS)) + ((INDXGLOB - 1L)%%NB) + 1L
  
  return( indx )
}

#' @export
g2lpair <- function(gi, gj, bldim, ICTXT)
{
  dum <- 0
  
  grid <- blacs(ICTXT=ICTXT)
  
  i <- indxg2l(gi, bldim[1L], 0L, 0L, grid$NPROW)
  j <- indxg2l(gj, bldim[2L], 0L, 0L, grid$NPCOL)
  
  return( c(i, j) )
}



#' @export
indxl2g <- function(INDXLOC, NB, IPROC, ISRCPROC, NPROCS)
{
  indx <- NPROCS*NB*(as.integer((INDXLOC - 1L)/NB)) + 
          (INDXLOC - 1L)%%NB + 
          ((NPROCS+IPROC-ISRCPROC)%%NPROCS)*NB + 
          1L
  
  return( indx )
}

#' @export
l2gpair <- function(i, j, bldim, ICTXT)
{
  grid <- blacs(ICTXT=ICTXT)
  
  gi <- indxl2g(i, bldim[1L], grid$MYROW, 0L, grid$NPROW)
  gj <- indxl2g(j, bldim[2L], grid$MYCOL, 0L, grid$NPCOL)
  
  return( c(gi, gj) )
}


