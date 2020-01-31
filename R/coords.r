#' Local to Global/Global to Local Indexing
#' 
#' Get the local index given global information.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param INDXGLOB
#' Global index.
#' @param INDXLOC
#' Local index.
#' @param NB
#' Block size.
#' @param IPROC
#' Coordinate of the process whose local info is to be determined.
#' @param ISRCPROC
#' The coordinate of the process that possesses the first row/column of the distributed matrix.  That's always 0 pbdDMAT.
#' @param NPROCS
#' Total number of processors over which matrix is distributed.
#' @return The local index.
#' 
#' @name coords
#' @rdname coords
#' @export
indxg2l <- function(INDXGLOB, NB, IPROC, ISRCPROC, NPROCS)
{
  indx <- NB*as.integer((INDXGLOB - 1L)/(NB*NPROCS)) + ((INDXGLOB - 1L)%%NB) + 1L
  indx
}

#' @name coords
#' @rdname coords
#' @export
indxl2g <- function(INDXLOC, NB, IPROC, ISRCPROC, NPROCS)
{
  indx <- NPROCS*NB*(as.integer((INDXLOC - 1L)/NB)) + 
          (INDXLOC - 1L)%%NB + 
          ((NPROCS+IPROC-ISRCPROC)%%NPROCS)*NB + 
          1L
  
  indx
}



#' Global to Local/Local to Global Pair Indexing
#' 
#' Get the local index-pair given global information.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param gi,gj
#' Global indices.
#' @param i,j
#' Local indices.
#' @param bldim
#' Blocking dimension
#' @param ICTXT
#' BLACS context.
#' @return The local index-pair.
#' 
#' @name coordspair
#' @rdname coordspair
#' @export
g2lpair <- function(gi, gj, bldim, ICTXT)
{
  dum <- 0
  
  grid <- blacs(ICTXT=ICTXT)
  
  i <- indxg2l(gi, bldim[1L], 0L, 0L, grid$NPROW)
  j <- indxg2l(gj, bldim[2L], 0L, 0L, grid$NPCOL)
  
  c(i, j)
}

#' @name coordspair
#' @rdname coordspair
#' @export
l2gpair <- function(i, j, bldim, ICTXT)
{
  grid <- blacs(ICTXT=ICTXT)
  
  gi <- indxl2g(i, bldim[1L], grid$MYROW, 0L, grid$NPROW)
  gj <- indxl2g(j, bldim[2L], grid$MYCOL, 0L, grid$NPCOL)
  
  c(gi, gj)
}



#' indxg2p
#' 
#' Computes the process coordinate which contains the entry of a
#' distributed matrix specified by a global index INDXGLOB.  
#' Simplified reimplementation of the ScaLAPACK aux INDXG2P function.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param INDXGLOB
#' Global index.
#' @param NB
#' Block size.
#' @param NPROCS
#' Total number of processors over which matrix is distributed.
#' @return The process coordinate.
#' 
#' @export
base.indxg2p <- function(INDXGLOB, NB, NPROCS)
{
  ISRCPROC <- 0L
  
  ret <- (ISRCPROC + as.integer((INDXGLOB - 1L) / NB)) %% NPROCS
  ret
}



#' numroc2
#' 
#' A better version of NUMROC (NUMber Rows Or Columns).  Returns the local
#' dimension given global matrix + distribution parameters.
#' 
#' For advanced users only. See pbdDMAT for high-level functions.
#' 
#' @param N
#' Global number of rows/cols.
#' @param NB
#' Block size.
#' @param IPROC
#' Coordinate of the process whose local info is to be determined.
#' @param NPROCS
#' Total number of processors over which matrix is distributed.
#' @return The local dimension.
#' 
#' @export
numroc2 <- function(N, NB, IPROC, NPROCS)
{
  ISRCPROC <- 0L
  
  MYDIST <- (NPROCS + IPROC - ISRCPROC) %% NPROCS
  NBLOCKS <- as.integer(N / NB)
  ldim <- as.integer(NBLOCKS / NPROCS) * NB
  EXTRABLKS <- NBLOCKS %% NPROCS
  
  if (is.na(EXTRABLKS))
    EXTRABLKS <- 0L
  
  if (MYDIST < EXTRABLKS)
    ldim <- ldim + NB
  else if (MYDIST == EXTRABLKS)
    ldim <- ldim + N %% NB
  
  ldim
}




#' Interchange Between Process Number and BLACS Coordinates
#' 
#' Grabs the existing BLACS context grid information.
#' 
#' For advanced users only. These functions are simple recreations of the BLACS
#' routines \code{BLACS_PNUM} and \code{BLACS_PCOORD}. The former gets the
#' process number associated with the BLACS process grid location
#' \code{c(MYPROW, MYPCOL)}, while the latter does the reverse.
#' 
#' @param ICTXT 
#' BLACS context number.
#' @param PROW,PCOL 
#' BLACS grid location row/column
#' @param PNUM
#' process rank
#' 
#' @return \code{pnum} returns an integer; \code{pcoord} returns a list
#' containing elements \code{PROW} and \code{PCOL}.
#' 
#' @keywords BLACS
#' 
#' @examples
#' spmd.code <- "
#'   suppressMessages(library(pbdMPI))
#'   suppressMessages(library(pbdBASE))
#'   init.grid()
#' 
#'   ### get the ICTXT = 0 BLACS coordsinates for process 3
#'   myCoords <- base.pcoord(ICTXT = 0, PNUM = 3)
#'   comm.print(myCoords)
#' 
#'   ### get the ICTXT = 1 BLACS coordsinates for process 3
#'   myCoords <- base.pcoord(ICTXT = 1, PNUM = 3)
#'   comm.print(myCoords)
#'
#'   ### get the ICTXT = 2 BLACS coordsinates for process 3
#'   myCoords <- base.pcoord(ICTXT = 2, PNUM = 3)
#'   comm.print(myCoords)
#' 
#'   finalize()
#' "
#' pbdMPI::execmpi(spmd.code = spmd.code, nranks = 4L)
#' 
#' @name pcoords
#' @rdname pcoords
#' @export
base.pnum <- function(ICTXT, PROW, PCOL)
{
  blacs_ <- base.blacs(ICTXT=ICTXT)
  NPCOL <- blacs_$NPCOL
  
  PNUM <- PROW * NPCOL + PCOL
  PNUM
}

pnum <- base.pnum



#' @rdname pcoords
#' @export
base.pcoord <- function(ICTXT, PNUM)
{
  blacs_ <- blacs(ICTXT=ICTXT)
  nprows <- blacs_$NPROW
  
  PROW <- as.integer(PNUM / nprows)
  PCOL <- PNUM %% nprows
  
  list(PROW=PROW, PCOL=PCOL)
}

pcoord <- base.pcoord



#' g2lcoord
#' 
#' Global to local coordinates with explicit ownership given.
#' 
#' @param dim
#' Global dimension.
#' @param bldim
#' Blocking dimension.
#' @param gi,gj
#' Global row and column indices, respectively.
#' @param gridinfo
#' The return of \code{base.blacs(ICTXT(x))}.  See the Details section
#' for more information.
#' 
#' @return
#' For the process that owns the desired local data at global indices
#' \code{(gi, gj)}, the return is the local index.  Otherwise, \code{NA}
#' is returned.
#' 
#' @useDynLib pbdBASE R_g2lcoord
#' @export
g2lcoord <- function(dim, bldim, gi, gj, gridinfo)
{
  .Call(R_g2lcoord, as.integer(dim), as.integer(bldim), 
        as.integer(gi), as.integer(gj), gridinfo)
}
