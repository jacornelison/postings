*
* $Id: mnseti.F,v 1.1 2006/11/10 17:49:36 tlc Exp $
*
* $Log: mnseti.F,v $
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
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
#include "minuit/pilot.h"
      SUBROUTINE MNSETI(TIT)
#include "minuit/d506dp.inc"
CC       Called by user to set or change title of current task.
CC
#include "minuit/d506cm.inc"
      CHARACTER*(*) TIT
      CTITL = TIT
      RETURN
      END
