C $Header: /u/gcmpack/MITgcm/pkg/pp81/pp81_diagnostics_init.F,v 1.1 2011/10/26 09:29:36 mlosch Exp $
C $Name:  $

#include "WALL_MODEL_OPTIONS.h"

CBOP
C     !ROUTINE: WALL_MODEL_DIAGNOSTICS_INIT
C     !INTERFACE:
      SUBROUTINE WALL_MODEL_DIAGNOSTICS_INIT( myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE WALL_MODEL_DIAGNOSTICS_INIT
C     | o Routine to initialize list of all available diagnostics
C     |   for WALL_MODEL package
C     *==========================================================*
C     \ev
C     !USES:
      IMPLICIT NONE

C     === Global variables ===

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_DIAGNOSTICS
C     !LOCAL VARIABLES :
C     === Local variables ===
C     diagNum   :: diagnostics number in the (long) list of available diag.
C     diagName  :: local short name (8c) of a diagnostics
C     diagCode  :: local parser field with characteristics of the diagnostics
C              cf head of S/R DIAGNOSTICS_INIT_EARLY or DIAGNOSTICS_MAIN_INIT
C     diagUnits :: local string (16c): physical units of a diagnostic field
C     diagTitle :: local string (80c): description of field in diagnostic
      INTEGER       diagNum
      CHARACTER*8   diagName
      CHARACTER*16  diagCode
      CHARACTER*16  diagUnits
      CHARACTER*(80) diagTitle

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

        diagName  = 'WMviscUr'
        diagTitle = 'WALL_MODEL Upt vertical eddy viscosity coefficient'
        diagUnits = 'm^2/s           '
        diagCode  = 'SU P    M1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'WMviscVr'
        diagTitle = 'WALL_MODEL Vpt vertical eddy viscosity coefficient'
        diagUnits = 'm^2/s           '
        diagCode  = 'SV P    M1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'WMdiffTS'
        diagTitle = 'WALL_MODEL Combined vertical diffusion coefficient
     &               for tracers'
        diagUnits = 'm^2/s           '
        diagCode  = 'SM P    M1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'WMdiffT '
        diagTitle = 'WALL_MODEL Temp vertical diffusion coefficient
     &               for tracers'
        diagUnits = 'm^2/s           '
        diagCode  = 'SM P    M1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'WMdiffS '
        diagTitle = 'WALL_MODEL Salt vertical diffusion coefficient
     &               for tracers'
        diagUnits = 'm^2/s           '
        diagCode  = 'SM P    M1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'WMCdU  '
        diagTitle = 'WALL_MODEL quadratic drag coefficient
     &               on U points' 
        diagUnits = 'm^2/s           '
        diagCode  = 'SU P    M1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'WMCdV  '
        diagTitle = 'WALL_MODEL quadratic drag coefficient
     &               on V points' 
        diagUnits = 'm^2/s           '
        diagCode  = 'SV P    M1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'WMCdC  '
        diagTitle = 'WALL_MODEL quadratic drag coefficient
     &               on C points' 
        diagUnits = 'm^2/s           '
        diagCode  = 'SM P    M1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'WMuModU '
        diagTitle = 'WALL_MODEL velocity on U points' 
        diagUnits = 'm^2/s           '
        diagCode  = 'SU      M1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'WMuModV '
        diagTitle = 'WALL_MODEL velocity on V points' 
        diagUnits = 'm^2/s           '
        diagCode  = 'SV      L1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'WMuModC '
        diagTitle = 'WALL_MODEL velocity on C points' 
        diagUnits = 'm^2/s           '
        diagCode  = 'SM      M1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'WMKEsurf'
        diagTitle = 'WALL_MODEL KE' 
        diagUnits = 'm^2/s           '
        diagCode  = 'SM      M1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I       diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END
