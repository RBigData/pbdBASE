!     PDLANGE subroutine wrapper for use with F77_CALL in C
      SUBROUTINE SUBPDLANGE(VALUE, NORM, M, N, A, IA, JA, DESCA, WORK)
      ! #######################################################################
      CHARACTER           NORM
      INTEGER             M, N, IA, JA, DESCA( 9 )
      DOUBLE PRECISION    VALUE, A( * ), WORK( * )
      DOUBLE PRECISION            PDLANGE
      
      VALUE = PDLANGE(NORM, M, N, A, IA, JA, DESCA, WORK)
      
      RETURN
      END
