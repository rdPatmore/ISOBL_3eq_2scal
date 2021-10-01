C $Header: /u/gcmpack/MITgcm/pkg/pp81/pp81_calc_visc.F,v 1.3 2009/10/08 20:08:20 jmc Exp $
C $Name:  $
#include "WALL_MODEL_OPTIONS.h"

      SUBROUTINE WALL_MODEL_CALC(
     I        bi,bj,myThid) 

CBOP
C     *==========================================================*
C     | SUBROUTINE WALL_MODEL_CALC                          |
C     | o Add contrubution to net viscosity from                 | 
C     |   WALL_MODEL mixing                                      |
C     *==========================================================*
      IMPLICIT NONE

C     == GLobal variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "WALL_MODEL.h"
#include "SHELFICE.h"
#include "DYNVARS.h"

C     == Routine arguments ==
C     bi, bj :: tile numbers
      INTEGER bi,bj
      INTEGER myThid

#ifdef ALLOW_WALL_MODEL

C     == Local variables ==
C     i, j, k    :: Loop counters
      INTEGER i,j,k
      _RL  uFld(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  vFld(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  KE(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  uSq, vSq
CEOP

      DO k=1,Nr

C     Make local copies of horizontal flow field
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
         uFld(i,j) = uVel(i,j,k,bi,bj)
         vFld(i,j) = vVel(i,j,k,bi,bj)
        ENDDO
       ENDDO

C     ---------------
C     O  U points
C     ---------------

       CALL WALL_MODEL_CALC_KE(bi,bj,k,selectKEscheme,kTopC,
     &                  uFld,vFld,KE,myThid)
       IF ( WALLMODELselectDragQuadr.EQ.0 ) THEN
C-     average grid-cell-center KE to get velocity norm @ U.pt
        DO j=1-OLy,sNy+OLy-1
         DO i=1-OLx+1,sNx+OLx-1
           uSq = 0. _d 0
           IF ( k.EQ.MAX( kTopC(i-1,j,bi,bj),kTopC(i,j,bi,bj) ) ) THEN
            uSq = KE(i,j)+KE(i-1,j)
            uModU(i,j,bi,bj) = SQRT( MAX(1. _d -6, uSq) )
           ENDIF
         ENDDO
        ENDDO
       ELSEIF ( WALLMODELselectDragQuadr.EQ.1 ) THEN
C-     calculate locally velocity norm @ U.pt (local U & 4 V averaged)
        DO j=1-OLy,sNy+OLy-1
         DO i=1-OLx+1,sNx+OLx-1
           uSq = 0. _d 0
           IF ( k.EQ.MAX( kTopC(i-1,j,bi,bj),kTopC(i,j,bi,bj) ) ) THEN
            uSq = uFld(i,j)*uFld(i,j)
     &        + ( (vFld(i-1, j )*vFld(i-1, j )*hFacS(i-1, j ,k,bi,bj)
     &            +vFld( i , j )*vFld( i , j )*hFacS( i , j ,k,bi,bj))
     &          + (vFld(i-1,j+1)*vFld(i-1,j+1)*hFacS(i-1,j+1,k,bi,bj)
     &            +vFld( i ,j+1)*vFld( i ,j+1)*hFacS( i ,j+1,k,bi,bj))
     &          )*recip_hFacW(i,j,k,bi,bj)*0.25 _d 0
            uModU(i,j,bi,bj) = SQRT( MAX(1. _d -6, uSq) )
           ENDIF
         ENDDO
        ENDDO
       ELSEIF ( WALLMODELselectDragQuadr.EQ.2 ) THEN
C-     same as above but using wet-point method to average 4 V
        DO j=1-OLy,sNy+OLy-1
         DO i=1-OLx+1,sNx+OLx-1
          uSq = 0. _d 0
c         IF (k .EQ. kTopC(i,j,bi,bj)) THEN
          IF ( k.EQ.MAX( kTopC(i-1,j,bi,bj),kTopC(i,j,bi,bj) ) ) THEN
           uSq = ( hFacS(i-1, j ,k,bi,bj) + hFacS( i , j ,k,bi,bj) )
     &         + ( hFacS(i-1,j+1,k,bi,bj) + hFacS( i ,j+1,k,bi,bj) )
           IF ( uSq.GT.zeroRL ) THEN
            uSq = uFld(i,j)*uFld(i,j)
     &        +( (vFld(i-1, j )*vFld(i-1, j )*hFacS(i-1, j ,k,bi,bj)
     &           +vFld( i , j )*vFld( i , j )*hFacS( i , j ,k,bi,bj))
     &         + (vFld(i-1,j+1)*vFld(i-1,j+1)*hFacS(i-1,j+1,k,bi,bj)
     &           +vFld( i ,j+1)*vFld( i ,j+1)*hFacS( i ,j+1,k,bi,bj))
     &        )/uSq
           ELSE
            uSq = uFld(i,j)*uFld(i,j)
           ENDIF
           uModU(i,j,bi,bj) = SQRT( MAX(1. _d -6, uSq) )
          ENDIF
         ENDDO
        ENDDO
       ENDIF

       DO j=1-OLy+1,sNy+OLy-1
        DO i=1-OLx,sNx+OLx-1
         IF ( k.EQ.MAX( kTopC(i-1,j,bi,bj),kTopC(i,j,bi,bj) ) ) THEN
          IF ( variableCd ) THEN
           BLDragCoeffU(i,j,bi,bj) = ( 1 / (recip_shiKarman *
     &                               LOG(0.5 * hFacW(i,j,k,bi,bj) *
     &                               drF(k) / roughness_length)) ) ** 2
          ELSE
           BLDragCoeffU(i,j,bi,bj) = WALLMODELDragQuadratic
          ENDIF
         ENDIF
        ENDDO
       ENDDO

C     ---------------
C     O  V Points
C     ---------------

       IF ( WALLMODELselectDragQuadr.EQ.0 ) THEN
C-     average grid-cell-center KE to get velocity norm @ U.pt
        DO j=1-OLy+1,sNy+OLy-1
         DO i=1-OLx,sNx+OLx-1
          vSq = 0. _d 0
          IF ( k.EQ.MAX( kTopC(i,j-1,bi,bj),kTopC(i,j,bi,bj) ) ) THEN
           vSq = KE(i,j)+KE(i,j-1)
           uModV(i,j,bi,bj) = SQRT( MAX(1. _d -6, vSq) )
          ENDIF
         ENDDO
        ENDDO
       ELSEIF ( WALLMODELselectDragQuadr.EQ.1 ) THEN
C-     calculate locally velocity norm @ U.pt (local U & 4 V averaged)
        DO j=1-OLy+1,sNy+OLy-1
         DO i=1-OLx,sNx+OLx-1
          vSq = 0. _d 0
          IF ( k.EQ.MAX( kTopC(i,j-1,bi,bj),kTopC(i,j,bi,bj) ) ) THEN
           vSq = vFld(i,j)*vFld(i,j)
     &       + ( (uFld( i ,j-1)*uFld( i ,j-1)*hFacW( i ,j-1,k,bi,bj)
     &           +uFld( i , j )*uFld( i , j )*hFacW( i , j ,k,bi,bj))
     &         + (uFld(i+1,j-1)*uFld(i+1,j-1)*hFacW(i+1,j-1,k,bi,bj)
     &           +uFld(i+1, j )*uFld(i+1, j )*hFacW(i+1, j ,k,bi,bj))
     &         )*recip_hFacS(i,j,k,bi,bj)*0.25 _d 0
           uModV(i,j,bi,bj) = SQRT( MAX(1. _d -6, vSq) )
          ENDIF
         ENDDO
        ENDDO
       ELSEIF ( WALLMODELselectDragQuadr.EQ.2 ) THEN
C-     same as above but using wet-point method to average 4 V
        DO j=1-OLy+1,sNy+OLy-1
         DO i=1-OLx,sNx+OLx-1
          vSq = 0. _d 0
          IF ( k.EQ.MAX( kTopC(i,j-1,bi,bj),kTopC(i,j,bi,bj) ) ) THEN
           vSq = ( hFacW( i ,j-1,k,bi,bj) + hFacW( i , j ,k,bi,bj) )
     &         + ( hFacW(i+1,j-1,k,bi,bj) + hFacW(i+1, j ,k,bi,bj) )
           IF ( vSq.GT.zeroRL ) THEN
            vSq = vFld(i,j)*vFld(i,j)
     &        +( (uFld( i ,j-1)*uFld( i ,j-1)*hFacW( i ,j-1,k,bi,bj)
     &           +uFld( i , j )*uFld( i , j )*hFacW( i , j ,k,bi,bj))
     &         + (uFld(i+1,j-1)*uFld(i+1,j-1)*hFacW(i+1,j-1,k,bi,bj)
     &           +uFld(i+1, j )*uFld(i+1, j )*hFacW(i+1, j ,k,bi,bj))
     &         )/vSq
           ELSE
            vSq = vFld(i,j)*vFld(i,j)
           ENDIF
           uModV(i,j,bi,bj) = SQRT( MAX(1. _d -6, vSq) )
          ENDIF
         ENDDO
        ENDDO
       ENDIF
       DO j=1-OLy+1,sNy+OLy-1
        DO i=1-OLx,sNx+OLx-1
         IF ( k.EQ.MAX( kTopC(i,j-1,bi,bj),kTopC(i,j,bi,bj) ) ) THEN
          IF ( variableCd ) THEN
C          The below formualtion is for a hydrolically smooth base
           BLDragCoeffV(i,j,bi,bj) = ( 1 / (recip_shiKarman *
     &                               LOG(0.5 * hFacS(i,j,k,bi,bj) *
     &                               drF(k) / roughness_length)) ) ** 2 
          ELSE
           BLDragCoeffV(i,j,bi,bj) = WALLMODELDragQuadratic
          ENDIF
         ENDIF
        ENDDO
       ENDDO

C     ---------------
C     O  C Points
C     ---------------

       DO j=1-OLy+1,sNy+OLy-1
        DO i=1-OLx,sNx+OLx-1
         IF (k .EQ. kTopC(i,j,bi,bj)) THEN
          IF ( variableCd ) THEN
           BLDragCoeffC(i,j,bi,bj) = ( 1 / (recip_shiKarman *
     &                                LOG(0.5 * hFacC(i,j,k,bi,bj) *
     &                                drF(k) / roughness_length)) ) ** 2
          ELSE
           BLDragCoeffC(i,j,bi,bj) = WALLMODELDragQuadratic
          ENDIF
          uModC(i,j,bi,bj) = SQRT( MAX(1. _d -6, 2 * KE(i,j)) ) 
         ENDIF
        ENDDO
       ENDDO
      ENDDO

#ifdef ALLOW_DIAGNOSTICS
       IF ( useDiagnostics ) THEN
         CALL DIAGNOSTICS_FILL(uModU,'WMuModU ',0,1,2,bi,bj,myThid)
         CALL DIAGNOSTICS_FILL(uModV,'WMuModV ',0,1,2,bi,bj,myThid)
         CALL DIAGNOSTICS_FILL(uModC,'WMuModC ',0,1,2,bi,bj,myThid)
         CALL DIAGNOSTICS_FILL(BLDragCoeffU,
     &                         'WMCdU   ',0,1,2,bi,bj,myThid)
         CALL DIAGNOSTICS_FILL(BLDragCoeffV,
     &                         'WMCdV   ',0,1,2,bi,bj,myThid)
         CALL DIAGNOSTICS_FILL(BLDragCoeffC,
     &                         'WMCdC   ',0,1,2,bi,bj,myThid)
       ENDIF
#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_WALL_MODEL */

      RETURN
      END
