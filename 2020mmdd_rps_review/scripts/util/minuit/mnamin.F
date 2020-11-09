*
* $Id: mnamin.F,v 1.1 2006/11/10 17:49:36 tlc Exp $
*
* $Log: mnamin.F,v $
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
      SUBROUTINE MNAMIN(FCN,FUTIL)
#include "minuit/d506dp.inc"
CC        Called  from many places.  Initializes the value of AMIN by
CC        calling the user function. Prints out the function value and
CC        parameter values if Print Flag value is high enough.
CC
#include "minuit/d506cm.inc"
      EXTERNAL FCN,FUTIL
      NPARX = NPAR
      IF (ISW(5) .GE. 1) WRITE (ISYSWR,'(/A,A)') ' FIRST CALL TO ',
     + 'USER FUNCTION AT NEW START POINT, WITH IFLAG=4.'
      CALL MNEXIN(X)
      CALL FCN(NPARX,GIN,FNEW,U,4,FUTIL)
      NFCN = NFCN + 1
      AMIN = FNEW
      EDM = BIGEDM
      RETURN
      END
