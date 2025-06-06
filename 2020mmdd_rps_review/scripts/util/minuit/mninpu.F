*
* $Id: mninpu.F,v 1.1 2006/11/10 17:49:36 tlc Exp $
*
* $Log: mninpu.F,v $
* Revision 1.1  2006/11/10 17:49:36  tlc
* added modified matmin and minuit subdir
*
* Revision 1.1  2005/11/08 16:27:41  pryke
* added stand alone minuit
*
* Revision 1.1  2005/11/08 16:14:06  pryke
* added stand alone minuit code
*
* Revision 1.1.1.1  2000/06/08 11:19:20  andras
* import of MINUIT from CERNlib 2000
*
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
#include "minuit/pilot.h"
      SUBROUTINE MNINPU(IUNIT,IERR)
#include "minuit/d506dp.inc"
CC      called by the user to SET INPUT to IUNIT,
CC      an alternative to MNSTIN where the user can specify just
CC      a logical unit number and he is not interrogated about
CC      open files and rewinding, all that is the responsibility
CC      of the user and cannot be fixed interactively.
#include "minuit/d506cm.inc"
C
      IERR = 0
C                              IUNIT = 0, revert to previous input file
      IF (IUNIT .EQ. 0) THEN
        IF (NSTKRD .EQ. 0)  THEN
           WRITE (ISYSWR, '(A)') ' CALL TO MNINPU(0) IGNORED'
           WRITE (ISYSWR, '(A)') ' ALREADY READING FROM PRIMARY INPUT'
        ELSE
          ISYSRD = ISTKRD(NSTKRD)
          NSTKRD = NSTKRD - 1
        ENDIF
C
C                               new input file
      ELSE
          IF (NSTKRD .GE. MAXSTK)  THEN
          WRITE (ISYSWR, '(A)') ' INPUT FILE STACK SIZE EXCEEDED.'
          GO TO 800
          ENDIF
        NSTKRD = NSTKRD + 1
        ISTKRD(NSTKRD) = ISYSRD
        ISYSRD = IUNIT
      ENDIF
C
      RETURN
  800 IERR = 1
      RETURN
      END
