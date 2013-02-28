! Copyright 2013, Schmidt


! Initialize BLACS communicator, process grid
      SUBROUTINE mpi_blacs_initialize(NPROW, NPCOL, ICTXT, MYROW, MYCOL)
      ! Inputs
      INTEGER NPROW, NPCOL
      ! Outputs
      INTEGER ICTXT, MYROW, MYCOL
      ! External
      EXTERNAL SL_INIT, BLACS_GRIDINFO
      
      CALL SL_INIT( ICTXT, NPROW, NPCOL )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

      RETURN
      END

