*
* $Id: mnbins.F,v 1.1 2006/11/10 17:49:36 tlc Exp $
*
* $Log: mnbins.F,v $
* Revision 1.1  2006/11/10 17:49:36  tlc
* added modified matmin and minuit subdir
*
* Revision 1.1  2005/11/08 16:27:41  pryke
* added stand alone minuit
*
* Revision 1.1  2005/11/08 16:14:06  pryke
* added stand alone minuit code
*
* Revision 1.1.1.1  2000/06/08 11:19:18  andras
* import of MINUIT from CERNlib 2000
*
* Revision 1.1.1.1  1996/03/07 14:31:28  mclareni
* Minuit
*
*
#include "minuit/pilot.h"
      SUBROUTINE MNBINS(A1,A2,NAA,BL,BH,NB,BWID)
#include "minuit/d506dp.inc"
C         SUBROUTINE TO DETERMINE REASONABLE HISTOGRAM INTERVALS
C         GIVEN ABSOLUTE UPPER AND LOWER BOUNDS  A1 AND A2
C         AND DESIRED MAXIMUM NUMBER OF BINS NAA
C         PROGRAM MAKES REASONABLE BINNING FROM BL TO BH OF WIDTH BWID
C         F. JAMES,   AUGUST, 1974 , stolen for Minuit, 1988
      PARAMETER (ZERO=0.0, ONE=1.0)
      AL = MIN(A1,A2)
      AH = MAX(A1,A2)
      IF (AL.EQ.AH)  AH = AL + 1.
C         IF NAA .EQ. -1 , PROGRAM USES BWID INPUT FROM CALLING ROUTINE
      IF (NAA .EQ. -1)  GO TO 150
   10 NA = NAA - 1
      IF (NA .LT. 1)  NA = 1
C          GET NOMINAL BIN WIDTH IN EXPON FORM
   20 AWID = (AH-AL)/FLOAT(NA)
      LOG = INT(DLOG10(DBLE(AWID)))
      IF (AWID .LE. ONE)  LOG=LOG-1
      SIGFIG = AWID * (10.00 **(-LOG))
C         ROUND MANTISSA UP TO 2, 2.5, 5, OR 10
      IF(SIGFIG .GT. 2.0)  GO TO 40
      SIGRND = 2.0
      GO TO 100
   40 IF (SIGFIG .GT. 2.5)  GO TO 50
      SIGRND = 2.5
      GO TO 100
   50 IF(SIGFIG .GT. 5.0)  GO TO 60
      SIGRND =5.0
      GO TO 100
   60 SIGRND = 1.0
      LOG = LOG + 1
  100 CONTINUE
      BWID = SIGRND*10.0**LOG
      GO TO 200
C         GET NEW BOUNDS FROM NEW WIDTH BWID
  150 IF (BWID .LE. ZERO)  GO TO 10
  200 CONTINUE
      ALB = AL/BWID
      LWID=ALB
      IF (ALB .LT. ZERO)  LWID=LWID-1
      BL = BWID*FLOAT(LWID)
      ALB = AH/BWID + 1.0
      KWID = ALB
      IF (ALB .LT. ZERO)  KWID=KWID-1
      BH = BWID*FLOAT(KWID)
      NB = KWID-LWID
      IF (NAA .GT. 5)  GO TO 240
      IF (NAA .EQ. -1)  RETURN
C          REQUEST FOR ONE BIN IS DIFFICULT CASE
      IF (NAA .GT. 1 .OR. NB .EQ. 1)  RETURN
      BWID =  BWID*2.0
       NB  = 1
       RETURN
  240 IF (2*NB .NE. NAA)  RETURN
      NA = NA + 1
      GO TO 20
      END
