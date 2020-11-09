*
* $Id: d506dp.inc,v 1.1.2.1 2009/09/11 16:59:21 pryke Exp $
*
* $Log: d506dp.inc,v $
* Revision 1.1.2.1  2009/09/11 16:59:21  pryke
* added missing minuit/minuit subdir
*
* Revision 1.1  2008/04/17 23:53:07  pryke
* *** empty log message ***
*
* Revision 1.1  2005/11/08 16:27:41  pryke
* added stand alone minuit
*
* Revision 1.1  2005/11/08 16:14:06  pryke
* added stand alone minuit code
*
* Revision 1.1.1.1  2000/06/08 11:19:21  andras
* import of MINUIT from CERNlib 2000
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
#ifndef CERNLIB_MINUIT_D506DP_INC
#define CERNLIB_MINUIT_D506DP_INC
*
*
* d506dp.inc
*
#if defined(CERNLIB_DOUBLE)
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#endif
#if !defined(CERNLIB_DOUBLE)
C ************ SINGLE PRECISION VERSION *************
#endif

#endif
