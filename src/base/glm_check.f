! This Source Code Form is subject to the terms of the Mozilla Public
! License, v. 2.0. If a copy of the MPL was not distributed with this
! file, You can obtain one at http://mozilla.org/MPL/2.0/.

! Copyright 2013, Schmidt


! Check the FAMILY and LINK arguments for valid/supported possibilities
    ! 0: no problem
    ! -1 FAMILY is invalid or unsupported
    ! -2 LINK is invalid or unsupported
      INTEGER FUNCTION GLM_CHECK_FAM_LINK(FAMILY, LINK)
      ! IN/OUT
      CHARACTER*8         FAMILY, LINK
      ! Parameters
      INTEGER             BAD_FAM, BAD_LINK
      PARAMETER ( BAD_FAM = -1, BAD_LINK = -2 )
      
      
      IF (FAMILY.EQ.'BINOMIAL') THEN
          IF (LINK.NE.'CLOGLOG'  .AND. 
     $        LINK.NE.'LOG'      .AND.
     $        LINK.NE.'LOGIT')    THEN
                  GLM_CHECK_FAM_LINK = BAD_LINK
          END IF
      
      ELSE IF (FAMILY.EQ.'GAMMA') THEN
          IF (LINK.NE.'IDENTITY'   .AND.
     $        LINK.NE.'LOG'        .AND. 
     $        LINK.NE.'INVERSE')   THEN
                  GLM_CHECK_FAM_LINK = BAD_LINK
          END IF
      
      ELSE IF (FAMILY.EQ.'GAUSSIAN') THEN
          IF (LINK.NE.'IDENTITY'      .AND.
     $        LINK.NE.'LOG'           .AND. 
     $        LINK.NE.'INVERSE')       THEN
                  GLM_CHECK_FAM_LINK = BAD_LINK
          END IF
      
      ELSE IF (FAMILY.EQ.'POISSON') THEN
          IF (LINK.NE.'IDENTITY'     .AND.
     $        LINK.NE.'LOG'          .AND. 
     $        LINK.NE.'SQRT')         THEN
                  GLM_CHECK_FAM_LINK = BAD_LINK
          END IF
      
      ELSE
          ! Family not supported
          CHECK_FAM_LINK = BAD_FAM
      END IF
      
      RETURN
      END

