C $Header: /u/gcmpack/MITgcm/pkg/pp81/pp81_calc_diff.F,v 1.6 2015/02/23 21:20:15 jmc Exp $
C $Name:  $

#include "WALL_MODEL_OPTIONS.h"

      SUBROUTINE WALL_MODEL_CALC_DIFF(
     I        bi,bj,iMin,iMax,jMin,jMax,
     I        wallGamma,
     U        KappaRx,KappaDiag,
     I        myThid)

C     *==========================================================*
C     | SUBROUTINE WALL_MODEL_CALC_DIFF                                |
C     | o Add contrubution to net diffusivity from WALL_MODEL mixing   |
C     *==========================================================*
      IMPLICIT NONE

C     == GLobal variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SHELFICE.h"
#ifdef ALLOW_3D_DIFFKR
# include "DYNVARS.h"
#endif
#ifdef ALLOW_WALL_MODEL
# include "WALL_MODEL.h"
#endif

C     == Routine arguments ==
C     bi, bj,   :: tile indices
C     iMin,iMax :: Range of points for which calculation is done
C     jMin,jMax :: Range of points for which calculation is done
C     KappaRx   :: vertical diffusivity array
C     myThid    :: Instance number for this innvocation of WALL_MODEL_CALC_DIFF

      INTEGER bi,bj,iMin,iMax,jMin,jMax
      _RL KappaRx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL KappaDiag(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL wallGamma
      INTEGER myThid

#ifdef ALLOW_WALL_MODEL

C     == Local variables ==
C     i,j,k     :: Loop counters
      INTEGER i,j,k

C--   Add Vertical diffusivity contribution from WALL_MODEL
C-    do all levels :
      DO j=jMin,jMax
       DO i=iMin,iMax
        k = kTopC(i,j,bi,bj)
        KappaRx(i,j,bi,bj) = SQRT(BLDragCoeffC(i,j,bi,bj)) *
     &                        uModC(i,j,bi,bj) * wallGamma
#ifdef ALLOW_DIAGNOSTICS
        IF ( useDiagnostics ) THEN
         KappaDiag(i,j) = KappaRx(i,j,bi,bj)
     &                  * hFacC(i,j,k,bi,bj) * drF(k) 
        ENDIF
#endif /* ALLOW_DIAGNOSTICS */
       ENDDO
      ENDDO
#endif /* ALLOW_WALL_MODEL */

      RETURN
      END
