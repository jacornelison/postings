*
* $Id: mninex.F,v 1.1 2006/11/10 17:49:36 tlc Exp $
*
* $Log: mninex.F,v $
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
      SUBROUTINE MNINEX(PINT)
#include "minuit/d506dp.inc"
CC        Transforms from internal coordinates (PINT) to external
CC        parameters (U).   The minimizing routines which work in
CC        internal coordinates call this routine before calling FCN.
#include "minuit/d506cm.inc"
      DIMENSION PINT(*)
      DO 100 J= 1, NPAR
      I = NEXOFI(J)
      IF (NVARL(I) .EQ. 1) THEN
         U(I) = PINT(J)
      ELSE
         U(I) = ALIM(I) + 0.5*(SIN(PINT(J)) +1.0) * (BLIM(I)-ALIM(I))
      ENDIF
  100 CONTINUE
      RETURN
      END
