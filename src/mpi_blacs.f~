! Initialize BLACS communicator, process grid
      SUBROUTINE mpi_blacs_initialize(NPROW, NPCOL, ICTXT, MYROW, MYCOL)
      EXTERNAL   SL_INIT, BLACS_GRIDINFO
      
      ! Inputs
      INTEGER NPROW, NPCOL
      ! Outputs
      INTEGER ICTXT, MYROW, MYCOL
      
      CALL SL_INIT( ICTXT, NPROW, NPCOL )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

      RETURN
      END

! Row/column-process sums
      SUBROUTINE row_col_sums( ICTXT, SCOPE, M, N, LDA, A )
      EXTERNAL   DGSUM2D
      
      ! Constants
      INTEGER           DEST
      PARAMETER       ( DEST = -1 )
      CHARACTER*1       TOP
      PARAMETER       ( TOP = ' ' )
      ! Inputs
      INTEGER           ICTXT, M, N, LDA
      CHARACTER*1       SCOPE
      ! In/Output
      DOUBLE PRECISION  A(M)
      
      CALL DGSUM2D( ICTXT, SCOPE, TOP, 
     $              M, N, A, LDA, DEST, DEST )

      RETURN
      END
