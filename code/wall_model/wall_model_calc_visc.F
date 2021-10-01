C $Header: /u/gcmpack/MITgcm/pkg/pp81/pp81_calc_visc.F,v 1.3 2009/10/08 20:08:20 jmc Exp $
C $Name:  $
#include "WALL_MODEL_OPTIONS.h"

      SUBROUTINE WALL_MODEL_CALC_VISC(
     I        bi,bj,iMin,iMax,jMin,jMax,k,
     U        KappaRU,KappaRV,
     I        myThid)

CBOP
C     *==========================================================*
C     | SUBROUTINE WALL_MODEL_CALC_VISC                                |
C     | o Add contrubution to net viscosity from WALL_MODEL mixing     |
C     *==========================================================*
      IMPLICIT NONE

C     == GLobal variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "GRID.h"
#include "SHELFICE.h"
#include "WALL_MODEL.h"

C     == Routine arguments ==
C     bi, bj :: tile numbers
C     iMin, iMax, jMin, jMax :: Range of points for which calculation
C     maskUp :: land/water mask for Wvel points (above tracer level)
C     myThid :: Instance number for this innvocation of WALL_MODEL_CALC_VISC
C
      INTEGER bi,bj,iMin,iMax,jMin,jMax,k
      _RL KappaRU(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      _RL KappaRV(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      INTEGER myThid

#ifdef ALLOW_WALL_MODEL

C     == Local variables ==
C     i, j    :: Loop counters
      INTEGER i,j
CEOP
C     This formualtion requires bottomVisc_pCell = .TRUE.
C     Discrepancy remains even with this flag.
C     Drag uses drC and this term uses drF.
      DO j=jMin,jMax
       DO i=iMin,iMax
        IF ( k.EQ.MAX( kTopC(i-1,j,bi,bj),kTopC(i,j,bi,bj) ) ) THEN
         KappaRU(i,j,k) = uModU(i,j,bi,bj)
     &        * BLDragCoeffU(i,j,bi,bj) * 
     &        * (0.5 * hFacW(i,j,k,bi,bj) * drF(k)) 
         nuU(i,j) = KappaRU(i,j,k)
#ifdef ALLOW_DIAGNOSTICS
         IF ( useDiagnostics ) THEN
          CALL DIAGNOSTICS_FILL(nuU,'WMviscUr',0,1,2,bi,bj,myThid)
         ENDIF
#endif /* ALLOW_DIAGNOSTICS */
        ENDIF
        IF ( k.EQ.MAX( kTopC(i,j-1,bi,bj),kTopC(i,j,bi,bj) ) ) THEN
         KappaRV(i,j,k) =  uModV(i,j,bi,bj)
     &       * BLDragCoeffV(i,j,bi,bj) * 
     &       * (0.5 * hFacS(i,j,k,bi,bj) * drF(k)) 
         nuV(i,j) = KappaRV(i,j,k)
#ifdef ALLOW_DIAGNOSTICS
         IF ( useDiagnostics ) THEN
          CALL DIAGNOSTICS_FILL(nuV,'WMviscVr',0,1,2,bi,bj,myThid)
         ENDIF
#endif /* ALLOW_DIAGNOSTICS */
        ENDIF
       ENDDO
      ENDDO

#endif /* ALLOW_WALL_MODEL */

      RETURN
      END
