*
* $Id: mnintr.F,v 1.1 2006/11/10 17:49:36 tlc Exp $
*
* $Log: mnintr.F,v $
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
      SUBROUTINE MNINTR(FCN,FUTIL)
#include "minuit/d506dp.inc"
CC       Called by user. Interfaces to MNREAD to allow user to change
CC       easily from Fortran-callable to interactive mode.
CC
#include "minuit/d506cm.inc"
      EXTERNAL FCN,FUTIL
      IFLGIN = 3
      CALL MNREAD(FCN,IFLGIN,IFLGUT,FUTIL)
      WRITE (ISYSWR,'(2A/)')  ' END OF MINUIT COMMAND INPUT. ',
     +      '   RETURN TO USER PROGRAM.'
      RETURN
      END
