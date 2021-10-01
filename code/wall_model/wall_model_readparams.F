C $Header: /u/gcmpack/MITgcm/pkg/pp81/pp81_readparms.F,v 1.4 2014/05/27 23:41:29 jmc Exp $
C $Name:  $

#include "WALL_MODEL_OPTIONS.h"

CBOP
C !ROUTINE: WALL_MODEL_READPARMS

C !INTERFACE: ==========================================================
      SUBROUTINE WALL_MODEL_READPARMS( myThid )

C !DESCRIPTION:
C     Initialize WALL_MODEL parameters, read in data.pp81

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "WALL_MODEL.h"
#include "SHELFICE.h"
#include "PARAMS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_WALL_MODEL

C !LOCAL VARIABLES: ====================================================
C  iUnit                :: unit number for I/O
C  msgBuf               :: message buffer
      INTEGER i,j,k
      INTEGER iUnit
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

      NAMELIST /WALL_MODEL_PARM01/
     &     recip_shiKarman,
     &     variableCd,
     &     roughness_length,
     &     WALLMODELselectDragQuadr,
     &     WALLMODELDragQuadratic

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      IF ( .NOT.useWallModel ) THEN
C-    pkg WALL_MODEL is not used
        _BEGIN_MASTER(myThid)
C-    Track pkg activation status:
         WallModelisON = .FALSE.
C     print a (weak) warning if data.wall_model is found
         CALL PACKAGES_UNUSED_MSG( 'useWallModel', ' ', ' ' )
        _END_MASTER(myThid)
        RETURN
      ENDIF

C This routine has been called by the main model so we set our
C internal flag to indicate we are in business
      WallModelisON = .TRUE.

C Set defaults values for parameters in WALL_MODEL.h
      recip_shiKarman  = 1. _d 0 / 0.4 _d 0
      variableCd       = .FALSE.
      roughness_length = 1. _d -2 / 30
      WALLMODELDragQuadratic = SHELFICEDragQuadratic
      WALLMODELSelectDragQuadr = UNSET_RL
      DO k=1,Nr+1
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
c         nuU(i,j,k) = viscArNr(k)
c         nuV(i,j,k) = viscArNr(k)
         nuU(i,j) = 1. _d -5
         nuV(i,j) = 1. _d -5
        ENDDO
       ENDDO
      ENDDO
   
C set the drag coefficient scheme
      IF (WALLMODELselectDragQuadr .EQ. UNSET_RL) THEN
       IF (SHELFICEselectDragQuadr .EQ. -1) THEN
        WALLMODELselectDragQuadr = 2
       ELSE
        WALLMODELselectDragQuadr = SHELFICEselectDragQuadr
       ENDIF
      ENDIF

C Open and read the data.wall_model file
      _BEGIN_MASTER(myThid)
      WRITE(msgBuf,'(A)') ' WALL_MODEL_READPARMS: opening 
     &                      data.wall_model'
      CALL PRINT_MESSAGE(msgBuf, standardMessageUnit,
     &                   SQUEEZE_RIGHT , 1)
      CALL OPEN_COPY_DATA_FILE(
     I                   'data.wall_model', 'WALL_MODEL_READPARMS',
     O                   iUnit,
     I                   myThid )
      READ(UNIT=iUnit,NML=WALL_MODEL_PARM01)
      WRITE(msgBuf,'(A)')
     &  ' WALL_MODEL_READPARMS: finished reading data.wall_model'
      CALL PRINT_MESSAGE(msgBuf, standardMessageUnit,
     &                   SQUEEZE_RIGHT , 1)

C Close the open data file
      CLOSE(iUnit)
      _END_MASTER(myThid)

C Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_WALL_MODEL */

      RETURN
      END
