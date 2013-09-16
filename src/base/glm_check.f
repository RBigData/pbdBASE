
! 0 is ok, -1 means family not supported, -2 means link not supported
      INTEGER FUNCTION CHECK_FAM_LINK(FAMILY, LINK)
      ! IN/OUT
      CHARACTER*8         FAMILY, LINK
      ! Parameters
      INTEGER             BAD_FAM, BAD_LINK
      PARAMETER ( BAD_FAM = -1, BAD_LINK = -2 )
      
      
      IF (FAMILY.EQ.'BINOMIAL') THEN
          IF (LINK.NE.'CLOGLOG'  .AND. 
     $        LINK.NE.'LOG'      .AND.
     $        LINK.NE.'LOGIT')    THEN
                  CHECK_FAM_LINK = BAD_LINK
          END IF
      
      ELSE IF (FAMILY.EQ.'GAMMA') THEN
          IF (LINK.NE.'IDENTITY'   .AND.
     $        LINK.NE.'LOG'        .AND. 
     $        LINK.NE.'INVERSE')   THEN
                  CHECK_FAM_LINK = BAD_LINK
          END IF
      
      ELSE IF (FAMILY.EQ.'GAUSSIAN') THEN
          IF (LINK.NE.'IDENTITY'      .AND.
     $        LINK.NE.'LOG'           .AND. 
     $        LINK.NE.'INVERSE')       THEN
                  CHECK_FAM_LINK = BAD_LINK
          END IF
      
      ELSE IF (FAMILY.EQ.'POISSON') THEN
          IF (LINK.NE.'IDENTITY'     .AND.
     $        LINK.NE.'LOG'          .AND. 
     $        LINK.NE.'SQRT')         THEN
                  CHECK_FAM_LINK = BAD_LINK
          END IF
      
      ELSE
          ! Family not supported
          CHECK_FAM_LINK = BAD_FAM
      END IF
      
      RETURN
      END




