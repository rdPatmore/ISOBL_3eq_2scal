C $Header: /u/gcmpack/MITgcm/pkg/shelfice/shelfice_update_masks.F,v 1.6 2014/09/11 19:20:38 jmc Exp $
C $Name: checkpoint65k $

#include "SHELFICE_OPTIONS.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif

CBOP
C     !ROUTINE: SHELFICE_UPDATE_MASKS
C     !INTERFACE:
      SUBROUTINE SHELFICE_UPDATE_MASKS(
     I     rF, recip_drF,
     U     hFacC,
     I     myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SHELFICE_UPDATE_MASKS
C     | o modify topography factor hFacC according to ice shelf
C     |   topography
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_SHELFICE
# include "SHELFICE.h"
#endif /* ALLOW_SHELFICE */

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     rF        :: R-coordinate of face of cell (units of r).
C     recip_drF :: Recipricol of cell face separation along Z axis ( units of r ).
C     hFacC     :: Fraction of cell in vertical which is open (see GRID.h)
C     myThid    :: Number of this instance of SHELFICE_UPDATE_MASKS
      _RS rF        (1:Nr+1)
      _RS recip_drF (1:Nr)
      _RS hFacC     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,1:Nr,nSx,nSy)
      INTEGER myThid

#ifdef ALLOW_SHELFICE
C     !LOCAL VARIABLES:
C     == Local variables ==
C     bi,bj   :: tile indices
C     I,J,K   :: Loop counters
      INTEGER bi, bj
      INTEGER I, J, K
      _RL hFacCtmp
      _RL hFacMnSz
CEOP

C     initialize R_shelfIce
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          R_shelfIce(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      IF ( SHELFICEtopoFile .NE. ' ' ) THEN
       _BARRIER
C Read the shelfIce draught using the mid-level I/O pacakage read_write_rec
C The 0 is the "iteration" argument. The 1 is the record number.
       CALL READ_REC_XY_RS( SHELFICEtopoFile, R_shelfIce,
     &      1, 0, myThid )
C Read the shelfIce draught using the mid-level I/O pacakage read_write_fld
C The 0 is the "iteration" argument. The ' ' is an empty suffix
C      CALL READ_FLD_XY_RS( SHELFICEtopoFile, ' ', R_shelfIce,
C    &      0, myThid )
      ENDIF
C- end setup R_shelfIce in the interior

C- fill in the overlap (+ BARRIER):
      _EXCH_XY_RS(R_shelfIce, myThid )

C     RDP 15/01/19 added to fix ice shelf problem
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO j=1,sNy
C        DO j=1-OLy,sNy+OLy
         IF (myPx .EQ. 1) THEN
          DO i=1-OLx,0
           R_shelfice(i,j,bi,bj) = R_shelfice(1,j,bi,bj)
          ENDDO
         ENDIF
         IF (myPx .EQ. nPx) THEN
          DO i=sNx+1,sNx+OLx
           R_shelfice(i,j,bi,bj) = R_shelfice(sNx,j,bi,bj)
          ENDDO
         ENDIF
        ENDDO
       ENDDO
      ENDDO

C--   Calculate lopping factor hFacC : Remove part outside of the domain
C     taking into account the Reference (=at rest) Surface Position Ro_shelfIce
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)

C--   compute contributions of shelf ice to looping factors
        DO K=1, Nr
         DO J=1-OLy,sNy+OLy
          DO I=1-OLx,sNx+OLx
C          RDP - Make w-e boundaries equal size in z when using p-cells
C          This is not right yet! under construction
           hFacMnSz=max( hFacMin, min(hFacMinDr*recip_drF(k),1. _d 0) )
           IF (myPx .EQ. 1) THEN
            IF (i .LT. OLx) THEN
             hFacMnSz = 1.0
            ENDIF
           ENDIF
           IF (myPx .EQ. nPx) THEN
            IF (i .GT. sNx-Olx) THEN
             hFacMnSz = 1.0
            ENDIF
           ENDIF
C      o Non-dimensional distance between grid boundary and model surface
           hFacCtmp = (rF(k)-R_shelfIce(I,J,bi,bj))*recip_drF(K)
C      o Reduce the previous fraction : substract the outside part.
           hFacCtmp = hFacC(I,J,K,bi,bj) - max( hFacCtmp, 0. _d 0)
C      o set to zero if empty Column :
           hFacCtmp = max( hFacCtmp, 0. _d 0)
C      o Impose minimum fraction and/or size (dimensional)
           IF (hFacCtmp.LT.hFacMnSz) THEN
            IF (hFacCtmp.LT.hFacMnSz*0.5) THEN
             hFacC(I,J,K,bi,bj)=0.
            ELSE
             hFacC(I,J,K,bi,bj)=hFacMnSz
            ENDIF
           ELSE
             hFacC(I,J,K,bi,bj)=hFacCtmp
           ENDIF
          ENDDO
         ENDDO
        ENDDO

#ifdef ALLOW_SHIFWFLX_CONTROL
C     maskSHI is a hack to play along with the general ctrl-package
C     infrastructure, where only the k=1 layer of a 3D mask is used
C     for 2D fields. We cannot use maskInC instead, because routines
C     like ctrl_get_gen and ctrl_set_unpack_xy require 3D masks.
        DO K=1,Nr
         DO J=1-OLy,sNy+OLy
          DO I=1-OLx,sNx+OLx
           maskSHI(I,J,K,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
        ENDDO
        DO K=1,Nr
         DO J=1-OLy,sNy+OLy
          DO I=1-OLx,sNx+OLx
           IF ( ABS(R_shelfice(I,J,bi,bj)) .GT. 0. _d 0
     &          .AND. hFacC(I,J,K,bi,bj) .NE. 0. _d 0 ) THEN
            maskSHI(I,J,K,bi,bj) = 1. _d 0
            maskSHI(I,J,1,bi,bj) = 1. _d 0
           ENDIF
          ENDDO
         ENDDO
        ENDDO
#endif /* ALLOW_SHIFWFLX_CONTROL */

C - end bi,bj loops.
       ENDDO
      ENDDO
#endif /* ALLOW_SHELFICE */

      RETURN
      END
