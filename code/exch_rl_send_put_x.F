C $Header: /u/gcmpack/MITgcm/eesupp/src/exch_rx_send_put_x.template,v 1.15 2012/09/06 15:25:01 jmc Exp $
C $Name: checkpoint65k $
#include "CPP_EEOPTIONS.h"

CBOP

C     !ROUTINE: EXCH_RL_SEND_PUT_X

C     !INTERFACE:
      SUBROUTINE EXCH_RL_SEND_PUT_X( array,
     I             myOLw, myOLe, myOLs, myOLn, myNz,
     I             exchWidthX, exchWidthY,
     I             thesimulationMode, thecornerMode, myThid )
      IMPLICIT NONE
C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE EXCH_RL_SEND_PUT_X
C     | o "Send" or "put" X edges for RL array.
C     *==========================================================*
C     | Routine that invokes actual message passing send or
C     | direct "put" of data to update X faces of an XY[R] array.
C     *==========================================================*

C     !USES:
C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "EESUPPORT.h"
#include "EXCH.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     array :: Array with edges to exchange.
C     myOLw :: West, East, North and South overlap region sizes.
C     myOLe
C     myOLn
C     myOLs
C     exchWidthX :: Width of data region exchanged.
C     exchWidthY
C     theSimulationMode :: Forward or reverse mode exchange ( provides
C                          support for adjoint integration of code. )
C     theCornerMode     :: Flag indicating whether corner updates are
C                          needed.
C     myThid            :: Thread number of this instance of S/R EXCH...
C     eBl               :: Edge buffer level
      INTEGER myOLw
      INTEGER myOLe
      INTEGER myOLs
      INTEGER myOLn
      INTEGER myNz
      _RL array(1-myOLw:sNx+myOLe,
     &          1-myOLs:sNy+myOLn,
     &          myNZ, nSx, nSy)
      INTEGER exchWidthX
      INTEGER exchWidthY
      INTEGER theSimulationMode
      INTEGER theCornerMode
      INTEGER myThid

C     !LOCAL VARIABLES:
C     == Local variables ==
C     i, j, k, iMin, iMax, iB    :: Loop counters and extents
C     bi, bj
C     biW, bjW                   :: West tile indices
C     biE, bjE                   :: East tile indices
C     eBl                        :: Current exchange buffer level
C     theProc, theTag, theType,  :: Variables used in message building
C     theSize
C     westCommMode               :: Working variables holding type
C     eastCommMode                  of communication a particular
C                                   tile face uses.
      INTEGER i, j, k, iMin, iMax, iB
      INTEGER bi, bj, biW, bjW, biE, bjE
      INTEGER eBl
      INTEGER westCommMode
      INTEGER eastCommMode


#ifdef ALLOW_USE_MPI
      INTEGER theProc, theTag, theType, theSize, mpiRc
# ifdef ALLOW_AUTODIFF_OPENAD_AMPI
      INTEGER mpiStatus(MPI_STATUS_SIZE)
      INTEGER pReqI
# endif
#endif
C--   Write data to exchange buffer
C     Various actions are possible depending on the communication mode
C     as follows:
C       Mode      Action
C     --------   ---------------------------
C     COMM_NONE  Do nothing
C
C     COMM_MSG   Message passing communication ( e.g. MPI )
C                Fill west send buffer from this tile.
C                Send data with tag identifying tile and direction.
C                Fill east send buffer from this tile.
C                Send data with tag identifying tile and direction.
C
C     COMM_PUT   "Put" communication ( UMP_, shmemput, etc... )
C                Fill east receive buffer of west-neighbor tile
C                Fill west receive buffer of east-neighbor tile
C                Sync. memory
C                Write data-ready Ack for east edge of west-neighbor
C                tile
C                Write data-ready Ack for west edge of east-neighbor
C                tile
C                Sync. memory
C
CEOP

     

#ifdef ALLOW_AUTODIFF_OPENAD_AMPI
# ifdef ALLOW_USE_MPI
      IF ( usingMPI ) THEN
      _BEGIN_MASTER(myThid)
       DO bj=1,nSy
        DO bi=1,nSx
          CALL ampi_awaitall (
     &         exchNReqsX(1,bi,bj) ,
     &         exchReqIdX(1,1,bi,bj) ,
     &         mpiStatus ,
     &         mpiRC )
        ENDDO
       ENDDO
      _END_MASTER(myThid)
      ENDIF
# endif
#endif

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C     Prevent anyone to access shared buffer while an other thread modifies it
      _BARRIER

C     Fill shared buffers from array values
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

        eBl = exchangeBufLevel(1,bi,bj)
        westCommMode = _tileCommModeW(bi,bj)
        eastCommMode = _tileCommModeE(bi,bj)
        biE = _tileBiE(bi,bj)
        bjE = _tileBjE(bi,bj)
        biW = _tileBiW(bi,bj)
        bjW = _tileBjW(bi,bj)

C     >>>>>>>>>>>>>>>>>>> FORWARD RUN <<<<<<<<<<<<<<<<<<<

        IF ( _theSimulationMode .EQ. FORWARD_SIMULATION ) THEN

C       o Send or Put west edge
         iMin = 1
         iMax = 1+exchWidthX-1
         IF ( westCommMode .EQ. COMM_MSG  ) THEN
C          PRINT *, 'RL west comm'
          iB = 0
          IF ( (verShift .NE. 0) .AND. (myNz .NE. 1) .AND.
     &       ( myPx .EQ. 1 ) ) THEN
           DO k=1 + verShift,myNz
            DO j=1,sNy
             DO i=iMin,iMax
              iB = iB + 1
              westSendBuf_RL(iB,eBl,bi,bj) = array(i,j,k,bi,bj)
             ENDDO
            ENDDO
           ENDDO
           DO k=myNz - verShift,myNz
            DO j=1,sNy
             DO i=iMin,iMax
              iB = iB + 1
             ENDDO
            ENDDO
           ENDDO
          ELSE
           DO k=1,myNz
            DO j=1,sNy
             DO i=iMin,iMax
              iB = iB + 1
              westSendBuf_RL(iB,eBl,bi,bj) = array(i,j,k,bi,bj)
             ENDDO
            ENDDO
           ENDDO
          ENDIF

         ELSEIF ( westCommMode .EQ. COMM_PUT  ) THEN
C          PRINT *, 'RL west put'
          iB  = 0
          IF ( (verShift .NE. 0) .AND. (myNz .NE. 1) .AND. 
     &         ( myPx .EQ. 1 ) ) THEN
C          DO k=1+verShift,myNz
           DO k=1 + verShift,myNz
            DO j=1,sNy
             DO i=iMin,iMax
              iB = iB + 1
              eastRecvBuf_RL(iB,eBl,biW,bjW) = array(i,j,k,bi,bj)
             ENDDO
            ENDDO
           ENDDO
           DO k=myNz - verShift,myNz
            DO j=1,sNy
             DO i=iMin,iMax
              iB = iB + 1
             ENDDO
            ENDDO
           ENDDO
          ELSE
           DO k=1,myNz
            DO j=1,sNy
             DO i=iMin,iMax
              iB = iB + 1
              eastRecvBuf_RL(iB,eBl,biW,bjW) = array(i,j,k,bi,bj)
             ENDDO
            ENDDO
           ENDDO
          ENDIF
         ELSEIF ( westCommMode .NE. COMM_NONE
     &    .AND.   westCommMode .NE. COMM_GET ) THEN
          STOP ' S/R EXCH: Invalid commW mode.'
         ENDIF

C       o Send or Put east edge
         iMin = sNx-exchWidthX+1
         iMax = sNx
         IF ( eastCommMode .EQ. COMM_MSG  ) THEN
C          PRINT *, 'RL east comm'
          iB = 0
          IF ( (verShift .NE. 0) .AND. (myNz .NE. 1) .AND. 
     &       ( myPx .EQ. nPx ) ) THEN
           DO k=1,verShift
            DO j=1,sNy
             DO i=iMin,iMax
              iB = iB + 1
             ENDDO
            ENDDO
           ENDDO
           DO k=1,myNz - verShift
            DO j=1,sNy
             DO i=iMin,iMax
              iB = iB + 1
              eastSendBuf_RL(iB,eBl,bi,bj) = array(i,j,k,bi,bj)
             ENDDO
            ENDDO
           ENDDO
          ELSE
           DO k=1,myNz
            DO j=1,sNy
             DO i=iMin,iMax
              iB = iB + 1
              eastSendBuf_RL(iB,eBl,bi,bj) = array(i,j,k,bi,bj)
             ENDDO
            ENDDO
           ENDDO
          ENDIF

         ELSEIF ( eastCommMode .EQ. COMM_PUT  ) THEN
C          PRINT *, 'RL east put'
          iB  = 0
          IF ( (verShift .NE. 0) .AND. (myNz .NE. 1) .AND.
     &       ( myPx .EQ. nPx ) ) THEN
           DO k=1,verShift
            DO j=1,sNy
             DO i=iMin,iMax
              iB = iB + 1
             ENDDO
            ENDDO
           ENDDO
           DO k=1,myNz - verShift
            DO j=1,sNy
             DO i=iMin,iMax
              iB = iB + 1
              westRecvBuf_RL(iB,eBl,biE,bjE) = array(i,j,k,bi,bj)
             ENDDO
            ENDDO
           ENDDO
          ELSE
           DO k=1,myNz
            DO j=1,sNy
             DO i=iMin,iMax
              iB = iB + 1
              westRecvBuf_RL(iB,eBl,biE,bjE) = array(i,j,k,bi,bj)
             ENDDO
            ENDDO
           ENDDO
          ENDIF
         ELSEIF ( eastCommMode .NE. COMM_NONE
     &    .AND.   eastCommMode .NE. COMM_GET  ) THEN
          STOP ' S/R EXCH: Invalid commE mode.'
         ENDIF

C     >>>>>>>>>>>>>>>>>>> ADJOINT RUN <<<<<<<<<<<<<<<<<<<
        ELSEIF ( _theSimulationMode .EQ. REVERSE_SIMULATION ) THEN

C       o Send or Put west edge
         iMin = 1-exchWidthX
         iMax = 0
         IF ( westCommMode .EQ. COMM_MSG  ) THEN
          iB = 0
          DO k=1,myNz
           DO j=1,sNy
            DO i=iMin,iMax
             iB = iB + 1
             westSendBuf_RL(iB,eBl,bi,bj) = array(i,j,k,bi,bj)
             array(i,j,k,bi,bj) = 0.0
            ENDDO
           ENDDO
          ENDDO
         ELSEIF ( westCommMode .EQ. COMM_PUT  ) THEN
          iB  = 0
          DO k=1,myNz
           DO j=1,sNy
            DO i=iMin,iMax
             iB = iB + 1
             eastRecvBuf_RL(iB,eBl,biW,bjW) = array(i,j,k,bi,bj)
             array(i,j,k,bi,bj) = 0.0
            ENDDO
           ENDDO
          ENDDO
         ELSEIF ( westCommMode .NE. COMM_NONE
     &    .AND.   westCommMode .NE. COMM_GET ) THEN
          STOP ' S/R EXCH: Invalid commW mode.'
         ENDIF

C       o Send or Put east edge
         iMin = sNx+1
         iMax = sNx+exchWidthX
         IF ( eastCommMode .EQ. COMM_MSG  ) THEN
          iB = 0
          DO k=1,myNz
           DO j=1,sNy
            DO i=iMin,iMax
             iB = iB + 1
             eastSendBuf_RL(iB,eBl,bi,bj) = array(i,j,k,bi,bj)
             array(i,j,k,bi,bj) = 0.0
            ENDDO
           ENDDO
          ENDDO
         ELSEIF ( eastCommMode .EQ. COMM_PUT  ) THEN
          iB  = 0
          DO k=1,myNz
           DO j=1,sNy
            DO i=iMin,iMax
             iB = iB + 1
             westRecvBuf_RL(iB,eBl,biE,bjE) = array(i,j,k,bi,bj)
             array(i,j,k,bi,bj) = 0.0
            ENDDO
           ENDDO
          ENDDO
         ELSEIF ( eastCommMode .NE. COMM_NONE
     &    .AND.   eastCommMode .NE. COMM_GET  ) THEN
          STOP ' S/R EXCH: Invalid commE mode.'
         ENDIF

        ENDIF

       ENDDO
      ENDDO

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
C--   Signal completetion ( making sure system-wide memory state is
C--                         consistent ).

C     ** NOTE ** We are relying on being able to produce strong-ordered
C     memory semantics here. In other words we assume that there is a
C     mechanism which can ensure that by the time the Ack is seen the
C     overlap region data that will be exchanged is up to date.
      IF ( exchNeedsMemSync  ) CALL MEMSYNC

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        eBl = exchangeBufLevel(1,bi,bj)
        biE = _tileBiE(bi,bj)
        bjE = _tileBjE(bi,bj)
        biW = _tileBiW(bi,bj)
        bjW = _tileBjW(bi,bj)
        westCommMode = _tileCommModeW(bi,bj)
        eastCommMode = _tileCommModeE(bi,bj)
        IF ( westCommMode.EQ.COMM_PUT ) eastRecvAck(eBl,biW,bjW) = 1
        IF ( eastCommMode.EQ.COMM_PUT ) westRecvAck(eBl,biE,bjE) = 1
        IF ( westCommMode.EQ.COMM_GET ) eastRecvAck(eBl,biW,bjW) = 1
        IF ( eastCommMode.EQ.COMM_GET ) westRecvAck(eBl,biE,bjE) = 1
       ENDDO
      ENDDO

C--   Make sure "ack" setting is seen system-wide.
C     Here strong-ordering is not an issue but we want to make
C     sure that processes that might spin on the above Ack settings
C     will see the setting.
C     ** NOTE ** On some machines we wont spin on the Ack setting
C     ( particularly the T90 ), instead we will use s system barrier.
C     On the T90 the system barrier is very fast and switches out the
C     thread while it waits. On most machines the system barrier
C     is much too slow and if we own the machine and have one thread
C     per process preemption is not a problem.
      IF ( exchNeedsMemSync  ) CALL MEMSYNC

C     Wait until all threads finish filling buffer
      _BARRIER

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#ifdef ALLOW_USE_MPI
      IF ( usingMPI ) THEN
C--   Send buffer data: Only Master Thread do proc communication
      _BEGIN_MASTER(myThid)

      DO bj=1,nSy
       DO bi=1,nSx

        eBl = exchangeBufLevel(1,bi,bj)
        westCommMode = _tileCommModeW(bi,bj)
        eastCommMode = _tileCommModeE(bi,bj)
        biE = _tileBiE(bi,bj)
        bjE = _tileBjE(bi,bj)
        biW = _tileBiW(bi,bj)
        bjW = _tileBjW(bi,bj)
        theType = _MPI_TYPE_RL
        theSize = sNy*exchWidthX*myNz

        IF ( westCommMode .EQ. COMM_MSG  ) THEN
C       Send buffer data (copied from west edge)
         theProc = tilePidW(bi,bj)
         theTag  = _tileTagSendW(bi,bj)
# ifndef ALLOW_AUTODIFF_OPENAD_AMPI
         exchNReqsX(1,bi,bj) = exchNReqsX(1,bi,bj)+1
         CALL MPI_Isend( westSendBuf_RL(1,eBl,bi,bj), theSize,
     &                   theType, theProc, theTag, MPI_COMM_MODEL,
     &                   exchReqIdX(exchNReqsX(1,bi,bj),1,bi,bj),
     &                   mpiRc )
# else
         pReqI=exchNReqsX(1,bi,bj)+1
         CALL ampi_isend_RL(
     &        westSendBuf_RL(1,eBl,bi,bj),
     &        theSize,
     &        theType,
     &        theProc,
     &        theTag,
     &        MPI_COMM_MODEL,
     &        exchReqIdX(pReqI,1,bi,bj),
     &        exchNReqsX(1,bi,bj),
     &        mpiStatus ,
     &        mpiRc )
# endif /* ALLOW_AUTODIFF_OPENAD_AMPI */
c        eastRecvAck(eBl,biW,bjW) = 1
        ENDIF

        IF ( eastCommMode .EQ. COMM_MSG  ) THEN
C       Send buffer data (copied from east edge)
         theProc = tilePidE(bi,bj)
         theTag  = _tileTagSendE(bi,bj)
# ifndef ALLOW_AUTODIFF_OPENAD_AMPI
         exchNReqsX(1,bi,bj) = exchNReqsX(1,bi,bj)+1
         CALL MPI_Isend( eastSendBuf_RL(1,eBl,bi,bj), theSize,
     &                   theType, theProc, theTag, MPI_COMM_MODEL,
     &                   exchReqIdX(exchNReqsX(1,bi,bj),1,bi,bj),
     &                   mpiRc )
# else
         pReqI=exchNReqsX(1,bi,bj)+1
         CALL ampi_isend_RL(
     &        eastSendBuf_RL(1,eBl,bi,bj) ,
     &        theSize ,
     &        theType ,
     &        theProc ,
     &        theTag ,
     &        MPI_COMM_MODEL ,
     &        exchReqIdX(pReqI,1,bi,bj) ,
     &        exchNReqsX(1,bi,bj),
     &        mpiStatus ,
     &        mpiRc )
# endif /* ALLOW_AUTODIFF_OPENAD_AMPI */
c        westRecvAck(eBl,biE,bjE) = 1
        ENDIF

       ENDDO
      ENDDO

      _END_MASTER(myThid)

      ENDIF
#endif /* ALLOW_USE_MPI */

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
      RETURN
      END
