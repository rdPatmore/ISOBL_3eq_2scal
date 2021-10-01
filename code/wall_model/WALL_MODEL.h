C $Header: /u/gcmpack/MITgcm/pkg/pp81/PP81.h,v 1.1 2004/09/02 09:11:54 mlosch Exp $
C $Name:  $

#ifdef ALLOW_WALL_MODEL

CBOP
C !ROUTINE: PP81.h

C !DESCRIPTION: \bv
C     /==========================================================\
C     | WALL_MODEL.h                                             |
C     | o Basic header for wall model of Jenkins (2019)          |
C     |   a 'law of the wall' parameterisation.                  |
C     |   Contains all wall_model field declarations.            |
C     \==========================================================/

C-----------------------------------------------------------------------
C
C Constants that can be set in data.pp
C     PPnRi - exponent of denominator of Eq(1) in PP81
C     PPviscMin, PPdiffMin   - minimum viscosity/diffusivity in
C                              surface layer.
C                              Only used if ALLOW_PP81MOD is defined
C     PPviscMax              - maximum allowed viscosity
C     PPnu0, PPalpha         - further parameters 
C     RiLimit                - minimum Richardson number, 
C                              follows from PPviscMax
C     PPdumpFreq, PPtaveFreq - analogue to dumpFreq and taveFreq
C                              (=default) 
C     PPmixingMaps - if true, include PP diagnostic maps in STDOUT
C     PPwriteState - if true, write PP state to file
C
C Time varying parameters computed by subroutine pp_calc
C     PPviscAr - Vertical eddy viscosity coefficient         (m^2/s)
C     PPdiffKr - Vertical diffusion coefficient for heat, 
C                salt and tracers                            (m^2/s)
C
C-----------------------------------------------------------------------
C \ev
CEOP

      _RL    recip_shiKarman, WALLMODELDragQuadratic
      _RL    WALLMODELselectDragQuadr, roughness_length
      COMMON /WALL_MODEL_PARMS_R/     
     &     recip_shiKarman, WALLMODELDragQuadratic,
     &     WALLMODElselectDragQuadr, roughness_length

      _RL uModU (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uModV (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL uModC (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL BLDragCoeffU (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL BLDragCoeffV (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL BLDragCoeffC (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL nuU(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL nuV(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      COMMON /WALL_MODEL_FIELDS/ 
     &     uModU, uModV, uModC,
     &     BLDragCoeffU, BLDragCoeffV, BLDragCoeffC,
     &     nuU, nuV

      LOGICAL wallModelisOn, variableCd
      COMMON /WALL_MODEL_PARMS_L/
     &     wallModelisOn, variableCd

#endif /* ALLOW_WALL_MODEL */
