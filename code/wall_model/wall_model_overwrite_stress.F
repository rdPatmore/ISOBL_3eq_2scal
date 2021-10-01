C $Header: /u/gcmpack/MITgcm/pkg/pp81/pp81_calc_visc.F,v 1.3 2009/10/08 20:08:20 jmc Exp $
C $Name:  $
#include "WALL_MODEL_OPTIONS.h"

      SUBROUTINE WALL_MODEL_OVERWRITE_STRESS(
     I        i,j,bi,bj,dirArg,
     U        dragCoeff)

CBOP
C     *==========================================================*
C     | SUBROUTINE WALL_MODEL_CALC_VISC                          |
C     | o Add contrubution to net viscosity from                 | 
C     |   WALL_MODEL mixing                                      |
C     *==========================================================*
      IMPLICIT NONE

C     == GLobal variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "GRID.h"
#include "WALL_MODEL.h"

C     == Routine arguments ==
C     bi, bj :: tile numbers
      INTEGER i,j,bi,bj
      INTEGER dirArg
      _RL dragCoeff(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

#ifdef ALLOW_WALL_MODEL

C     == Local variables ==
C     i, j, k    :: Loop counters
CEOP
      IF ( dirArg .EQ. 0) THEN  
       dragCoeff(i,j) = BLDragCoeffU(i,j,bi,bj)
      ELSEIF ( dirArg .EQ. 1) THEN
       dragCoeff(i,j) = BLDragCoeffV(i,j,bi,bj)
      ENDIF

#endif /* ALLOW_WALL_MODEL */

      RETURN
      END
