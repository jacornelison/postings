/****************************************************************************
This library should be self contained, and includes all tools needed to 
pass between horizon coords (az,el) and equatorial coords (ra,dec) and 
vice versa.  Effects included are: Precession,Nutation,Aberration.

ACCURACY: (measured against Xephem) is better than 0.2 arcseconds for source
positions.  Accuracy of Sun position (ra,dec) is good to about 10
arcseconds.  There is no altitude correction for az,el calculations, so the
accuracy of the Sun az,el will be affected - but it's a small effect.

here is a listing of the routines and their function

get_jd() : returns julian day given ctime
jd2lst() : returns Local Sidereal Time given Julian Day and Longitude
precess(): precesses ra,dec from eqnx1 to eqnx2.  works on vector inputs of 
           length 'npts'.
radec2azel(): converts RA,DEC to az,el given julian day,lat,lon.
              ra and dec are converted to az and el in place, resectively.
              RA,DEC inputs should be J2000 [DEGREES], and are precessed
              to the epoch implied by the julian day and lon.
              OUTPUT AZ IS SPECIFIED West of south (ie FSS frame, I think...)
azel2radec(): inverse of the above operation, WITH THE EXCEPTION OF THE 
              AZ DEFINITION. 180 should be added(subtracted) to(from)
              the output of radec2azel() to get the inverse.  Az here is
              specified E of North, as is standard (but not with boom).
nutate(): nutates.  called by *2*.
cor_nutate(): does the nutation correction.  called by *2*
cor_aberration(): does the aberration correction.  called by *2*
sunpos(): given julian day, calculates sun RA,DEC and elongation.
moonpos(): given julian day, calculates moon RA,DEC.
sixty(): given a double scalar, it returns the three arguments of
         the sexigesimal representation

last updated Jan 2004

-Bill Jones

Added coord.c from <http://dsnra.jpl.nasa.gov/pkg/astro-1.1/html/astro/>
hcc 2/4/2005

*****************************************************************************/
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <malloc.h>

#ifndef PI
#define PI 3.14159265358979323846264338327950288
void mmmult(double **,double **,int);

void mmmult(double **a,double **b,int n) {
  int i,j,k;
  double tmp[n][n];
  for (i=0;i<n;i++) for (j=0;j<n;j++) tmp[i][j]=0;

  for (i=0;i<n;i++) for (j=0;j<n;j++)
    for (k=0;k<n;k++) 
      tmp[i][j]+=((double)a[i][k]*(double)b[k][j]);
  for (i=0;i<n;i++) for (j=0;j<n;j++)
    a[i][j]=tmp[i][j];
}
#endif
#define PI_2 6.283185307

#define MOON_NCOEFF 60
#define MOON_MAXBUFF 1024

double get_jd(double t);
double jd2lst(double jd,double lon);
void precess(double *ra,double *dec,double eqnx1,double eqnx2,long npts);
void radec2azel(double *ra,double *dec,double *jd,double *lat,double *lon,long npts);
void azel2radec(double *az,double *el,double *jd,double *lat,double *lon,long npts);
void nutate(double jd,double *dpsi,double *deps);
void cor_nutate(double ra,double dec,double jd,double *dra,double *ddec,double *eps,double *dpsi);
void sunpos(double jd,double *ra,double *dec,double *sunelon);
void cor_aberration(double ra,double dec,double jd,double eps,double *dra,double *ddec);
double gcirc(double ra1,double ra2,double dec1,double dec2);
void sixty(double n,int *h,int *m,double *s);
void moonpos(double jd, double *raout, double *decout);
void mvmult(double **,double *,int);
void coord (double long_of_new_origin, double lat_of_new_origin,
	    double long_of_new_pole, double lat_of_new_pole,
            double old_long, double old_lat,
            double *new_long, double *new_lat);

// BEGIN SOURCE
// Code for calculating RA, DEC and back.  Good to about 0.3 arcsec, based on
// XEphem.

// stolen subroutine from wcjlib.h (hcc 1/28/2005)
void mvmult(double **a,double *b,int n) {
  int i,j;
  double tmp[n];
  for (i=0;i<n;i++) tmp[i]=0;
  for (i=0;i<n;i++) for (j=0;j<n;j++) tmp[i]+=a[i][j]*b[j];
  for (i=0;i<n;i++) b[i]=tmp[i];
}

// return julian day based on ctime
double get_jd(double t){
  return(2440587.5+t/(3600.0*24.0));
}

// calculates Local Sidereal Time from julian day, longitude [DEGREES]
// lst is returned in [hours]
// from Jean Meeus, "Astronomical Algorithms"

double jd2lst(double jd,double lon) {
  double c[]={280.46061837, 360.98564736629, 0.000387933, 38710000.0 };
  double jd2k=2451545.0;
  double t0=(jd-jd2k);
  double t=t0/36525.0;
  double gmt,lst;
  
  gmt=c[0] + (c[1] * t0) + t*t*(c[2] - t/ c[3]);
  //  lst=(gmt-lon)/15.0; //convert degrees to hours
  //  hcc, 10/17/2005
  lst=(gmt+lon)/15.0; //convert degrees to hours
  if (lst < 0.0) lst=24.+fmod(lst,24);
  else lst=fmod(lst,24);
  return lst;
}

// precess ra,dec from equinox 1 to equinox 2
// input RA and DEC are in DEGREES
void precess(double *rain,double *decin,double eq1,double eq2,long npts) {
  int i;
  double t,st,*ra,*dec,*ang,**m;
  double c0,c1,c2,s0,s1,s2;
  double d2r=(PI/180.0);
  double s2r=(d2r/3600.0);

  ra=(double *)malloc((size_t)npts*sizeof(double));
  dec=(double *)malloc((size_t)npts*sizeof(double));
  for (i=0;i<npts;i++) {
    ra[i]=d2r*rain[i];
    dec[i]=d2r*decin[i];
  }
  ang=(double *)malloc((size_t)3*sizeof(double));
  m=(double **)malloc((size_t)npts*sizeof(double));
  for (i=0;i<3;i++) m[i]=(double *)malloc((size_t)3*sizeof(double));

  t=0.001*(eq2-eq1);
  st=0.001*(eq1-2000.0);

  ang[0]=s2r*t*(23062.181 + st*(139.656 +0.0139*st) + t*(30.188 - 0.344*st+17.998*t));
  ang[1]=s2r*t*t*(79.280 + 0.410*st + 0.205*t) + ang[0];
  ang[2]=s2r*t*(20043.109 - st*(85.33 + 0.217*st)+t*(-42.665 - 0.217*st -41.833*t));
  c0=cos(ang[0]); s0=sin(ang[0]);
  c1=cos(ang[1]); s1=sin(ang[1]);
  c2=cos(ang[2]); s2=sin(ang[2]);

  m[0][0]=c0*c1*c2-s0*s1;
  m[0][1]=(-c0*s1-s0*c1*c2);
  m[0][2]=(-c1*s2);

  m[1][0]=s0*c1+c0*s1*c2;
  m[1][1]=c0*c1-s0*s1*c2;
  m[1][2]=(-s1*s2);

  m[2][0]=c0*s2;
  m[2][1]=(-s0*s2);
  m[2][2]=c2;

  for (i=0;i<npts;i++) {
    ang[0]=cos(dec[i])*cos(ra[i]);
    ang[1]=cos(dec[i])*sin(ra[i]);
    ang[2]=sin(dec[i]);
    mvmult(m,ang,3);
    rain[i]=(atan2(ang[1],ang[0])/d2r);
    decin[i]=(asin(ang[2])/d2r);
    while (rain[i]<0) rain[i]+=360.0;
  }

  free(ra);
  free(dec);
  free(ang);
  for (i=0;i<3;i++) free(m[i]);
  free(m);
 
}

//input RA,DEC J200 (FK5), julian day, and lat,lon [degrees]
//ra and dec [DEGREES] pointers replaced with az,el [DEGREES]
//note the az+=180 at the end - this is to match the convention
//used in the old source Az's, and my not be desireable in practice
//
void radec2azel(double *ra,double *dec,double *jd,double *lat,double *lon,long npts) {
  long i;
  double dtmp,ha,x,y,z,r;
  double d2r=(PI/180.0);
  double ch,sh,cd,sd,cl,sl;
  double dra,dra2,ddec,ddec2,eps,dpsi;

  for (i=0;i<npts;i++) {
    dtmp=(jd[i]-2451545)/365.25+2000;
    precess(&ra[i],&dec[i],2000,dtmp,1);

    cor_nutate(ra[i],dec[i],jd[i],&dra,&ddec,&eps,&dpsi);
    cor_aberration(ra[i],dec[i],jd[i],eps,&dra2,&ddec2);

    ra[i]+=((dra + dra2)/3600.);
    dec[i]+=((ddec + ddec2)/3600.);

    dtmp=jd2lst(jd[i],lon[i])*15.0;

    dtmp+=(dpsi*cos(eps)/3600.); //add correction in degrees

    if ((ha=dtmp-ra[i])<0) ha+=360.0;
    ha=fmod(ha,360.0);

    ch=cos(d2r*ha); sh=sin(d2r*ha);
    cd=cos(d2r*dec[i]); sd=sin(d2r*dec[i]);
    cl=cos(d2r*lat[i]); sl=sin(d2r*lat[i]);

    x = - ch * cd * sl + sd * cl ;
    y = - sh * cd;
    z = ch * cd * cl + sd * sl;
    r = sqrt(x*x + y*y);

    ra[i]=atan2(y,x)/d2r;
    dec[i]=atan2(z,r)/d2r;

    ra[i]+=180.0;
    while (ra[i] < -180.0) ra[i]+=360.;
    while (ra[i] > 180.0) ra[i]-=360.;
  }

}


//Note that this wants the az in the standard convention, ie not the output of
//radec2azel.  to be clear, to perform the inverse of radec2azel you 
// would want something like the following
//
//  radec2azel(&ra,&dec,&jd,&lat,&lon,1);
//  ra-=180.;
//  azel2radec(&ra,&dec,&jd,&lat,&lon,1);
//  
//  if ra and dec are initially J2000 (as they should be) this example
//  will result in returning RA and DEC as they were
//
void azel2radec(double *az,double *el,double *jd,double *lat,double *lon,long npts) {
  long i;
  double dtmp,ha;
  double d2r=(PI/180.0);
  double sa,se,sl,ca,ce,cl;
  double dra,dra2,ddec,ddec2,eps,dpsi;

  for (i=0;i<npts;i++) {
 
    cor_nutate(45.0,45.0,jd[i],&dra,&ddec,&eps,&dpsi);
    //just to get eps for ha correction
    sa=sin(d2r*az[i]); ca=cos(d2r*az[i]);
    se=sin(d2r*el[i]); ce=cos(d2r*el[i]);
    sl=sin(d2r*lat[i]); cl=cos(d2r*lat[i]);
 
    ha=atan2(-sa*ce,-ca*sl*ce+se*cl)/d2r;
    ha=fmod(ha,360.0);
    el[i]=asin(sl*se+cl*ce*ca)/d2r;
    
    dtmp=jd2lst(jd[i],lon[i])*15.0;
    dtmp+=(dpsi *cos(eps)/3600.);

    az[i]=fmod((dtmp-ha+360.0),360.0);
    
    cor_nutate(az[i],el[i],jd[i],&dra,&ddec,&eps,&dpsi);
    cor_aberration(az[i],el[i],jd[i],eps,&dra2,&ddec2);

    az[i]-=((dra+dra2)/3600.);
    el[i]-=((ddec+ddec2)/3600.);

    dtmp=(jd[i]-2451545)/365.25+2000;
    precess(&az[i],&el[i],dtmp,2000,1);
  }

}

void nutate(double jd,double *dpsi,double *deps) {
  int i,n=63;

  double c1[]={297.85036,445267.111480,-0.0019142,1.0/189474.0};
  double c2[]={357.52772,35999.0503400,-0.0001603,-1.0/300000.0};
  double c3[]={134.96298,477198.867398,0.0086972,1.0/56250.0};
  double c4[]={93.271910,483202.017538,-0.0036825,-1.0/327270.0};
  double c5[]={125.04452,-1934.1362610,0.00207080,1.0/450000.0};

  double d_lng[]={0,-2,0,0,0,0,-2,0,0,-2,-2,-2,0,2,0,2,0,0,-2,0,2,0,0,-2,0,-2,0,0,2,-2,0,-2,0,0,2,2,0,-2,0,2,2,-2,-2,2,2,0,-2,-2,0,-2,-2,0,-1,-2,1,0,0,-1,0,0,2,0,2};
  
  double m_lng[] ={0,0,0,0,1,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,2,1,0,-1,0,0,0,1,1,-1,0,0,0,0,0,0,-1,-1,0,0,0,1,0,0,1,0,0,0,-1,1,-1,-1,0,-1};

  double mp_lng[]={0,0,0,0,0,1,0,0,1,0,1,0,-1,0,1,-1,-1,1,2,-2,0,2,2,1,0,0,-1,0,-1,0,0,1,0,2,-1,1,0,1,0,0,1,2,1,-2,0,1,0,0,2,2,0,1,1,0,0,1,-2,1,1,1,-1,3,0};

  double f_lng[]={0,2,2,0,0,0,2,2,2,2,0,2,2,0,0,2,0,2,0,2,2,2,0,2,2,2,2,0,0,2,0,0,0,-2,2,2,2,0,2,2,0,2,2,0,0,0,2,0,2,0,2,-2,0,0,0,2,2,0,0,2,2,2,2};

  double om_lng[]={1,2,2,2,0,0,2,1,2,2,0,1,2,0,1,2,1,1,0,1,2,2,0,2,0,0,1,0,1,2,1,1,1,0,1,2,2,0,2,1,0,2,1,1,1,0,1,1,1,1,1,0,0,0,0,0,2,0,0,2,2,2,2};

  double sin_lng[]={-171996,-13187,-2274,2062,1426,712,-517,-386,-301,217,-158, 129, 123, 63, 63, -59, -58, -51, 48, 46, -38, -31, 29, 29, 26, -22,21, 17, 16, -16, -15, -13, -12, 11, -10, -8, 7, -7, -7, -7,6,6,6,-6,-6,5,-5,-5,-5,4,4,4,-4,-4,-4,3,-3,-3,-3,-3,-3,-3,-3};

 double sdelt[]={-174.2, -1.6, -0.2, 0.2, -3.4, 0.1, 1.2, -0.4, 0, -0.5, 0, 0.1,0,0,0.1, 0,-0.1,0,0,0,0,0,0,0,0,0,0, -0.1, 0, 0.1,  0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0}; 

  double cos_lng[]={92025, 5736, 977, -895, 54, -7, 224, 200, 129, -95,0,-70,-53,0, -33, 26, 32, 27, 0, -24, 16,13,0,-12,0,0,-10,0,-8,7,9,7,6,0,5,3,-3,0,3,3,0,-3,-3,3,3,0,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  double cdelt[] = {8.9, -3.1, -0.5, 0.5, -0.1, 0.0, -0.6, 0.0, -0.1, 0.3, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  double d,m,mp,f,omega;
  double t=(jd-2451545.0)/36525.0;
  double arg,sarg,carg;
  double tmp,d2r=(PI/180.0);

  tmp=0; d=0; m=0; mp=0; f=0; omega=0;  

  tmp=c1[0]+c1[1]*t+c1[2]*t*t+c1[3]*t*t*t;
  while (tmp < 0.0) tmp+=360.0;
  while (tmp > 360.0) tmp-=360.0;
  d=tmp*d2r;

  tmp=c2[0]+c2[1]*t+c2[2]*t*t+c2[3]*t*t*t;
  while (tmp < 0.0) tmp+=360.0;
  while (tmp > 360.0) tmp-=360.0;
  m=tmp*d2r;

  tmp=c3[0]+c3[1]*t+c3[2]*t*t+c3[3]*t*t*t;
  while (tmp < 0.0) tmp+=360.0;
  while (tmp > 360.0) tmp-=360.0;
  mp=tmp*d2r;

  tmp=c4[0]+c4[1]*t+c4[2]*t*t+c4[3]*t*t*t;
  while (tmp < 0.0) tmp+=360.0;
  while (tmp > 360.0) tmp-=360.0;
  f=tmp*d2r;

  tmp=c5[0]+c5[1]*t+c5[2]*t*t+c5[3]*t*t*t;
  while (tmp < 0.0) tmp+=360.0;
  while (tmp > 360.0) tmp-=360.0;
  omega=tmp*d2r;

  *dpsi=0; *deps=0;  

  for (i=0;i<n;i++) {
    arg=d_lng[i]*d + m_lng[i]*m +mp_lng[i]*mp + f_lng[i]*f +om_lng[i]*omega;
    sarg=sin(arg);
    carg=cos(arg);

    *dpsi+=(sdelt[i]*t+sin_lng[i])*sarg;
    *deps+=(cdelt[i]*t+cos_lng[i])*carg;
  }
  *dpsi*=0.0001; *deps*=0.0001;
}

//Input RA,DEC in DEGREES
void cor_nutate(double ra,double dec,double jd,double *dra,double *ddec,double *eps,double *dpsi) {
  
  double t=(jd-2451545.0)/36525.0;
  //time in Julian centuries from 2000.0
  double d_eps,eps0;
  double d2r=(PI/180.0);
  double d2as=(PI/(180.0*3600.0));
  double ce,se,x,y,z,x2,y2,z2,r,xyproj,ra2,dec2;

  nutate(jd,dpsi,&d_eps);

  eps0 = 23.4392911*3600.0 - 46.8150*t - 0.00059*t*t + 0.001813*t*t*t;
  *eps = (eps0 + d_eps)/3600.*d2r;
  ce = cos(*eps);
  se = sin(*eps);

  x = cos(ra*d2r) * cos(dec*d2r);
  y = sin(ra*d2r) * cos(dec*d2r);
  z = sin(dec*d2r);

  x2 = x - (y*ce + z*se)*(*dpsi)*d2as;
  y2 = y + (x*ce*(*dpsi) - z*d_eps)*d2as;
  z2 = z + (x*se*(*dpsi) + y*d_eps)*d2as;

  r = sqrt(x2*x2 + y2*y2 + z2*z2);
  xyproj = sqrt(x2*x2 + y2*y2);

  ra2=0.0; dec2=0.0;
 
  if ((xyproj==0.0) && (z!=0.0)) {
    //point at NCP/SCP
    dec2=asin(z2/r);
    ra2=0.0;
  }
  if (xyproj!=0) {
    //point is not at NCP/SCP
    ra2=atan2(y2,x2);
    dec2=asin(z2/r);
  }

  ra2/=d2r;
  dec2/=d2r;

  while (ra2<0) ra2+=360.0;

  // Return ra and dec corrections in arcseconds
  *dra = (ra2 - ra) * 3600.0;
  *ddec = (dec2 - dec) * 3600.0;

}

//SUNPOS: oblit and sunelon, ra and dec are left in DEGREES
void sunpos(double jd,double *ra,double *dec,double *sunelon) {
  double d2r=(PI/180.0);
  double t = (jd - 2415020.0)/36525.0; 
  //time in Julian centuries from 1900.0
  double l,me,mv,mm,mj,d;
  double ellcor,vencorr,marscorr,jupcorr,mooncorr;
  double longterm,omega,oblt;

  l=(279.696678+fmod((36000.768925*t),360.0))*3600.0;

  me= 358.475844 + fmod((35999.049750*t),360.0);
  ellcor= (6910.1 - 17.2*t)*sin(me*d2r) + 72.3*sin(2.0*me*d2r);
  l+=ellcor;

  //Venus perturbations using the mean anomaly of Venus MV

  mv = 212.603219 + fmod((58517.803875*t),360.0); 
  vencorr = 4.8 * cos((299.1017 + mv - me)*d2r) + 5.5 * cos((148.3133 +  2.0 * mv  -  2.0 * me )*d2r) +2.5 * cos((315.9433 +  2.0 * mv  -  3.0 * me )*d2r) +1.6 * cos((345.2533 +  3.0 * mv  -  4.0 * me )*d2r) + 1.0 * cos((318.15   +  3.0 * mv  -  5.0 * me )*d2r);
  l+= vencorr;

  //Mars perturbations using the mean anomaly of Mars MM

  mm = 319.529425 + fmod((19139.858500 * t),360.0);
  marscorr = 2.0 * cos((343.8883-2.0*mm+2.0*me)*d2r)+1.8*cos((200.4017-2.0*mm+me)*d2r);
  l+= marscorr;

  //Jupiter perturbations using the mean anomaly of Jupiter MJ

  mj = 225.328328 + fmod((3034.6920239*t),360.0);
  jupcorr = 7.2 * cos(( 179.5317 - mj + me )*d2r) +2.6 * cos((263.2167  -  mj)*d2r) +2.7 * cos(( 87.1450 -  2.0 * mj  +  2.0 * me ) *d2r) + 1.6 * cos((109.4933 -  2.0 * mj  +  me ) *d2r);
  l+= jupcorr;

  //Moons perturbations using the mean elongation of the Moon from Sun D

  d = 350.7376814  + fmod((445267.11422*t),360.0);
  mooncorr=6.5 * sin(d*d2r);
  l+=mooncorr;

  longterm = 6.4*sin(( 231.19+20.20*t)*d2r);
  l+= longterm;
  l=fmod(( l + 2592000.0),1296000.0); 
  *sunelon= l/3600.0;

  l-=20.5;

  // Nutation using the longitude of the Moons mean node OMEGA

  omega = 259.183275 - fmod(( 1934.142008* t ),360.0);
  l-=(17.2*sin(omega*d2r));

  oblt = 23.452294 - 0.0130125*t + (9.2*cos(omega*d2r))/3600.0;

  //Form Right Ascension and Declination
  l/=3600.0;
  *ra=(atan2(sin(l*d2r)*cos(oblt*d2r),cos(l*d2r))/d2r);
  while (*ra<0.0) *ra+=360;

  *dec=(asin(sin(l*d2r)*sin(oblt*d2r))/d2r);
}

void cor_aberration(double ra,double dec,double jd,double eps,double *dra,double *ddec) {
  double d2r =(PI/180.);
  double t = (jd -2451545.0)/36525.0; //julian centuries from J2000 of jd.
  double e,peri,k,sunra,sundec,sunlon;
  double cd,ce,cp,cs,ca;
  double sd,te,sp,ss,sa;
  double t1,t2,t3,t4;

  sunpos(jd,&sunra,&sundec,&sunlon);

  // Earth's orbital eccentricity
  e = 0.016708634 - 0.000042037*t - 0.0000001267*t*t;
  //longitude of perihelion, in degrees 
  peri = 102.93735 + 1.71946*t + 0.00046*t*t ;
  k=20.49552; //constant of aberration, in arcseconds

  cd = cos(dec*d2r);    sd = sin(dec*d2r);
  ce = cos(eps);        te = tan(eps);
  cp = cos(peri*d2r);   sp = sin(peri*d2r);
  cs = cos(sunlon*d2r); ss = sin(sunlon*d2r);
  ca = cos(ra*d2r);     sa = sin(ra*d2r);

  t1 = (ca*cs*ce+sa*ss)/cd;
  t2 = (ca*cp*ce+sa*sp)/cd;
  t3 = (cs*ce*(te*cd-sa*sd)+ca*sd*ss);
  t4 = (cp*ce*(te*cd-sa*sd)+ca*sd*sp);

  *dra = -k * t1 + e*k * t2;
  *ddec = -k * t3 + e*k * t4;
}

double gcirc(double ra1,double ra2,double dec1,double dec2) {
  double dis,radiff;

  radiff=fabs(ra1-ra2);
  if (radiff>PI) radiff=2*PI-radiff;
  dis=sin(dec1)*sin(dec2)+cos(dec1)*cos(dec2)*cos(radiff);
  if (fabs(dis) > 1) {
    fprintf(stderr,"GCIRC: DIS ERROR!\n"); 
    dis=0;
  } else dis=acos(dis);
  
  return dis;
}

void sixty(double n,int *h,int *m,double *s) {
  double tmp;
  if (n>0.0) {
    *h=(int)floor(n);
    *m=(int)floor(modf(n,&tmp)*60.0);
    *s=modf(modf(n,&tmp)*60,&tmp)*60.0;
  } else {
    n=fabs(n);
    *h=(-1)*(int)floor(n);
    *m=(int)floor(modf(n,&tmp)*60.0);
    *s=modf(modf(n,&tmp)*60,&tmp)*60.0;
  }
}

/*      FUNCTIONAL DESCRIPTION:  This subroutine converts the longitude-like
        and latitude-like coordinates of a point on a sphere
        into the corresponding coordinates in a different coordinate
        system that is specified by the coordinates of its origin
        and its north pole in the original coordinate system.
 
     REVISION HISTORY:

        DATE      PROGRAMMER      DESCRIPTION OF CHANGES 
        --------  --------------  -----------------------------------
         1Aug93   T.B.H.Kuiper    added explanation
        23Mar93   P.Withington    converted to C
        15Oct87   T.B.H.Kuiper    converted to COMAL
        24Nov82   T.B.H.Kuiper    HEAD documentation 
         6Jun79   T.B.H.Kuiper    last pre-HEAD revision 

      NOTES:  From Manchester and Gordon DOPSET

        Examples of use: 

        1. HOUR ANGLE, DECLINATION --> AZIMUTH, ELEVATION
           coord(PI, PI/2-lat, 0., lat, ha, dec, *az, *el)
           ----------------------------------------------------------
        2. AZIMUTH, ELEVATION --> HOUR ANGLE, DECLINATION
           coord(PI, PI/2-lat, 0., lat, az, el, *ha, *dec)
           ----------------------------------------------------------
        3. RIGHT ASCENSION, DECLINATION --> GALACTIC LONG & LAT
           ra_gal_pole   = (12. + 49./60.)*PI/12.
           dec_gal_pole  =  27.4          *PI/180.
           ra_gal_origin = (17. + 42.4/60)*PI/12.
           dec_gal_origin=-(28.0+ 55./60.)*PI/180.
                 see Blaauw, Gum, Pawsey and Westerhout,
                     MNRAS, 121, 123 (1960)
           coord(ra_gal_origin, dec_gal_origin,
                 ra_gal_pole, dec_gal_pole, ra, dec, *l, *b)
          --------------------------------------------------------------
        4. GALACTIC LONG & LAT --> RIGHT ASCENSION, DECLINATION
                 In general, whenever we know the forward transformation
                 (i.e. example 3 above) we may do the reverse transforma-
                 tion with at most two preliminary calls to coord to  cal-
                 culate the coordinates in system 2 of the pole and origin
                 which are given in system 1.
           coord(ra_gal_origin, dec_gal_origin,
                 ra_gal_pole, dec_gal_pole,
                 0., PI/2,
                 l_cel_pole, b_cel_pole)
           coord(ra_gal_origin, dec_gal_origin,
                 ra_gal_pole, dec_gal_pole,
                 0., 0.,
                 l_cel_origin, b_cel_origin)
           The results of this transformation are:
                 l,b of celestial pole   =  123.0013,   27.4000
                 l,b of celestial origin =   97.7435,  -60.1810
           Then:
           coord(l_cel_origin, b_cel_origin, l_cel_pole, b_cel_pole,
                 l, b, *ra, *dec)
    
        If long_of_new_pole, lat_of_new_pole, and long_of_new_origin are
        known, then there are only two values possible for lat_of_new_origin
        (except for degenerate special cases). The possible values of
        lat_of_new_origin may be calculated from

        sin(lat_of_new_origin)=
                  [sin(lat_of_new_pole)
                   +/- 2.*cos(lat_of_new_pole)^2
                         *cos(long_of_new_pole-long_of_new_origin)
                         *SQR{1.+cos(long_of_new_pole-long_of_new_origin)^2}]
                / [sin(lat_of_new_pole)^2
                     + {cos(lat_of_new_pole)
                        *cos(long_of_new_pole-long_of_new_origin)}^2]
    
        If long_of_new_pole, lat_of_new_pole, and lat_of_new_origin are known,
        the two possible values of long_of_new_origin may be calculated from

        cos(long_of_new_pole-long_of_new_origin)=
                               [1-sin(lat_of_new_origin)*sin(lat_of_new_pole)]
                                /[cos(lat_of_new_origin)*cos(lat_of_new_pole)]
    
 
        If, instead of long_of_new_origin and lat_of_new_origin, the longitude
        of the ascending node is known in both the old (an1) and new (an2)
        coordinate systems, then long_of_new_origin and lat_of_new_origin may
        be calculated by a preliminary call to coord:

        coord(0., 0., an1-long_of_new_pole, lat_of_new_pole,
              -an2, 0, *long_of_new_origin, *lat_of_new_origin)
   
*/
void coord (double long_of_new_origin, double lat_of_new_origin,
	    double long_of_new_pole, double lat_of_new_pole,
            double old_long, double old_lat,
            double *new_long, double *new_lat) {

/* long_of_new_origin;	       longitude of new origin in old system */
/* lat_of_new_origin;          latitude      "          "     "      */
/* long_of_new_pole;           longitude of new pole    "     "      */
/* lat_of_new_pole;            latitude      "          "     "      */
/* old_long;                   longitude of star        "     "      */
/* old_lat;                    latitude      "          "     "      */
/* *new_long;                  longitude     "          " new system */
/* *new_lat;                   latitude      "          "     "      */

	double sbp, cbp, sb0, sb1, cb1, sb2, cb2;
	double saa, caa, cb0, cbb, sbb, sa2, ca2;


	sb0 = sin (lat_of_new_origin);
	cb0 = cos (lat_of_new_origin);

	sbp = sin (lat_of_new_pole);
	cbp = cos (lat_of_new_pole);

	sb1 = sin (old_lat);
	cb1 = cos (old_lat);

	sb2 = sbp * sb1 + cbp * cb1 * cos (long_of_new_pole - old_long);

	*new_lat = asin (sb2);
	cb2 = cos (*new_lat);

	if ((cb2 != 0.0) && (cbp != 0.0))
	{
		saa = sin (long_of_new_pole - old_long) * cb1 / cb2;
		caa = (sb1 - sb2 * sbp) / (cb2 * cbp);
		cbb = sb0 / cbp;
		sbb = sin (long_of_new_pole - long_of_new_origin) * cb0;
		sa2 = (saa * cbb) - (caa * sbb);
		ca2 = (caa * cbb) + (saa * sbb);
		*new_long = atan2 (sa2, ca2);
		if (*new_long < 0.0)
			*new_long = *new_long + PI_2;
	}
	else
	{
		printf ("Error - bypassing AZ calculation. . .\n");
		*new_long = 0.;
	}
}

//Moon stuff

typedef struct pterm {
  int emod;
  long sincoeff;
  long coscoeff;
  int arg[4];
} PTerm;

PTerm slsr[MOON_MAXBUFF] = {
   { 0, 6288774, -20905355 , { 0, 0,  1, 0 },}, 
   { 0, 1274027,  -3699111 , { 2, 0, -1, 0 },}, 
   { 0,  658314,  -2955968 , { 2, 0,  0, 0 },}, 
   { 0,  213618,   -569925 , { 0, 0,  2, 0 },}, 
   { 0, -185116,     48888 , { 0, 1,  0, 0 },}, 
   { 0, -114332,     -3149 , { 0, 0,  0, 2 },}, 
   { 0,   58793,    246158 , { 2, 0, -2, 0 },}, 
   { 0,   57066,   -152138 , { 2, -1, -1, 0 },}, 
   { 0,   53322,   -170733 , { 2, 0,  1, 0 },}, 
   { 0,   45758,   -204586 , { 2, -1, 0, 0 },}, 
   { 0,  -40923,   -129620 , { 0,  1, -1, 0 },}, 
   { 0,  -34720,    108743 , { 1,  0, 0, 0 },}, 
   { 0,  -30383,    104755 , { 0,  1, 1, 0 },}, 
   { 0,   15327,     10321 , { 2,  0, 0, -2 },}, 
   { 0,  -12528,         0 , { 0,  0, 1, 2 },}, 
   { 0,   10980,     79661 , { 0,  0, 1, -2 },}, 
   { 0,   10675,    -34782 , { 4,  0, -1, 0 },}, 
   { 0,   10034,    -23210 , { 0,  0, 3, 0 },}, 
   { 0,    8548,    -21636 , { 4, 0, -2, 0 },}, 
   { 0,   -7888,     24208 , { 2, 1, -1, 0 },}, 
   { 0,   -6766,     30824 , { 2, 1,  0, 0 },}, 
   { 0,   -5163,     -8379 , { 1, 0, -1, 0 },}, 
   { 0,    4987,    -16675 , { 1, 1,  0, 0 },}, 
   { 0,    4036,    -12831 , { 2, -1, 1, 0 },}, 
   { 0,    3994,    -10445 , { 2, 0,  2, 0 },}, 
   { 0,    3861,    -11650 , { 4, 0,  0, 0 },}, 
   { 0,    3665,     14403 , { 2, 0, -3, 0 },}, 
   { 0,   -2689,     -7003 , { 0, 1, -2, 0 },}, 
   { 0,   -2602,         0 , { 2, 0, -1, 2 },}, 
   { 0,    2390,     10056 , { 2, -1, -2, 0 },}, 
   { 0,   -2348,      6322 , { 1, 0,  1, 0 },}, 
   { 0,    2236,     -9884 , { 2, -2, 0, 0 },}, 
   { 0,   -2120,      5751 , { 0, 1,  2, 0 },}, 
   { 0,   -2069,         0 , { 0, 2,  0, 0 },}, 
   { 0,    2048,     -4950 , { 2, -2, -1, 0 },}, 
   { 0,   -1773,      4130 , { 2, 0, 1, -2 },}, 
   { 0,   -1595,         0 , { 2, 0, 0, 2 },}, 
   { 0,    1215,     -3958 , { 4, -1, -1, 0 },}, 
   { 0,   -1110,         0 , { 0,  0,  2, 2 },}, 
   { 0,    -892,      3258 , { 3, 0, -1, 0 },}, 
   { 0,    -810,      2616 , { 2, 1,  1, 0 },}, 
   { 0,     759,     -1897 , { 4, -1, -2, 0 },}, 
   { 0,    -713,     -2117 , { 0,  2, -1, 0 },}, 
   { 0,    -700,      2354 , { 2, 2, -1, 0 },}, 
   { 0,     691,         0 , { 2, 1, -2, 0 },}, 
   { 0,     596,         0 , { 2, -1, 0, -2 },}, 
   { 0,     549,     -1423 , { 4, 0, 1, 0 },}, 
   { 0,     537,     -1117 , { 0, 0, 4, 0 },}, 
   { 0,     520,     -1571 , { 4, -1, 0, 0 },}, 
   { 0,    -487,     -1739 , { 1, 0, -2, 0 },}, 
   { 0,    -399,         0 , { 2, 1,  0, -2 },}, 
   { 0,    -381,     -4421 , { 0, 0, 2, -2 },}, 
   { 0,    -351,         0 , { 1, 1,  1, 0 },}, 
   { 0,    -340,         0 , { 3, 0, -2, 0 },}, 
   { 0,     330,         0 , { 4, 0, -3, 0 },}, 
   { 0,     327,         0 , { 2, -1, 2, 0 },}, 
   { 0,    -323,      1165 , { 0,  2, 1, 0 },}, 
   { 0,     299,         0 , { 1, 1, -1, 0 },}, 
   { 0,     294,         0 , { 2, 0,  3, 0 },}, 
   { 0,       0,      8752 , { 2, 0, -1, -2 }, }
  };

PTerm sb[MOON_MAXBUFF] = {
   { 0, 5128122 , 0, { 0, 0,  0,  1 },},  
   { 0,  280602 , 0, { 0, 0,  1,  1 },},  
   { 0,  277693 , 0, { 0,  0,  1, -1 },},  
   { 0,  173237 , 0, { 2,  0,  0, -1 },},  
   { 0,   55413 , 0, { 2,  0, -1,  1 },},  
   { 0,   46271 , 0, { 2,  0, -1, -1 },},  
   { 0,   32573 , 0, { 2,  0,  0,  1 },},  
   { 0,   17198 , 0, { 0,  0,  2,  1 },},  
   { 0,    9266 , 0, { 2,  0,  1, -1 },},  
   { 0,    8822 , 0, { 0,  0,  2,  1 },},  
   { 1,    8216 , 0, { 2, -1,  0, -1 },},  
   { 0,    4324 , 0, { 2,  0, -2, -1 },},  
   { 0,    4200 , 0, { 2,  0,  1,  1 },},  
   { 1,   -3359 , 0, { 2,  1,  0, -1 },},  
   { 1,    2463 , 0, { 2, -1, -1,  1 },},  
   { 1,    2211 , 0, { 2, -1,  0,  1 },},  
   { 1,    2065 , 0, { 2, -1, -1, -1 },},  
   { 1,   -1870 , 0, { 0,  1, -1, -1 },},  
   { 0,    1828 , 0, { 4,  0, -1, -1 },},  
   { 1,   -1794 , 0, { 0,  1,  0,  1 },},  
   { 0,   -1749 , 0, { 0,  0,  0,  3 },},  
   { 1,   -1565 , 0, { 0,  1, -1,  1 },},  
   { 0,   -1491 , 0, { 1,  0,  0,  1 },},  
   { 1,   -1475 , 0, { 0,  1,  1,  1 },},  
   { 1,   -1410 , 0, { 0,  1,  1, -1 },},  
   { 1,   -1344 , 0, { 0,  1,  0, -1 },},  
   { 0,   -1335 , 0, { 1,  0,  0, -1 },},  
   { 0,    1107 , 0, { 0,  0,  3,  1 },},  
   { 0,    1021 , 0, { 4,  0,  0, -1 },},  
   { 0,     833 , 0, { 4,  0, -1,  1 },},  
   { 0,     777 , 0, { 0,  0,  1, -3 },},  
   { 0,     671 , 0, { 4,  0, -2,  1 },},  
   { 0,     607 , 0, { 2,  0,  0, -3 },},  
   { 0,     596 , 0, { 2,  0,  2, -1 },},  
   { 1,     491 , 0, { 2, -1,  1, -1 },},  
   { 0,    -451 , 0, { 2,  0, -2,  1 },},  
   { 0,     439 , 0, { 0,  0,  3, -1 },},  
   { 0,     422 , 0, { 2,  0,  2,  1 },},  
   { 0,     421 , 0, { 2,  0, -3, -1 },},  
   { 1,    -366 , 0, { 2,  1, -1,  1 },},  
   { 1,    -351 , 0, { 2,  1,  0,  1 },},  
   { 0,     331 , 0, { 4,  0,  0,  1 },},  
   { 1,     315 , 0, { 2, -1,  1,  1 },},  
   { 2,     302 , 0, { 2, -2,  0, -1 },},  
   { 0,    -283 , 0, { 0,  0,  1,  3 },},  
   { 1,    -229 , 0, { 2,  1,  1, -1 },},  
   { 1,     223 , 0, { 1,  1,  0, -1 },},  
   { 1,     223 , 0, { 1,  1,  0, 1 },},  
   { 1,    -220 , 0, { 0,  1, -2, -1 },},  
   { 1,    -220 , 0, { 2,  1, -1, -1 },},  
   { 0,    -185 , 0, { 1,  0,  1,  1 },},  
   { 1,     181 , 0, { 2, -1, -2, -1 },},  
   { 1,    -177 , 0, { 0,  1,  2,  1 },},  
   { 0,     176 , 0, { 4,  0, -2, -1 },},  
   { 1,     166 , 0, { 4, -1, -1, -1 },},  
   { 0,    -164 , 0, { 1,  0,  1, -1 },},  
   { 0,     132 , 0, { 4,  0,  1, -1 },},  
   { 0,    -119 , 0, { 1,  0, -1, -1 },},  
   { 1,     115 , 0, { 4, -1,  0, -1 },},  
   { 2,     107 , 0, { 2, -2,  0,  1 },},  
};

typedef struct mparam
{
  double ra;
  double decl;

} OrbPar;

int getmoonpos(double T, OrbPar *MoonPar);

typedef struct sidtime
{
  double T;
  double JD;

} ATime;

double deg2rad(double ang) {
  double result;

  result = (ang * PI )/ 180.0;  
  return result;
}

double dot(double a[], double b[], int len) {
  double sum;
  int i;

  sum = 0.0;
  for ( i = 0; i < len ; i++ ) {
    sum += a[i]*b[i];
  }
  return(sum);
}

int getmoonpos(double T, OrbPar *MoonPar) {
  int i, j;
  int nterms;
  
  double ls, l, m, ms, d, f;
  double tmpsb[MOON_MAXBUFF], tmpsl[MOON_MAXBUFF], dmmf[MOON_MAXBUFF];
  
  double e, Om;
  double a1, a2, a3;
  double dotsb, dotslsr;

  double slsum, srsum, sbsum;
  double slres, srres, sbres;
  double dpsi, eps, eps0, deps;
  double lambda, beta, sigma, delta;
  double Nom, Den;
  double U;
  double ra, decl, sindecl; 
  

  //Updated as in second edition of Meeus..  wcj 2000

  // OLD VERSION
  /* mean length of moon */
  //ls = 218.3164591 + 481267.88134236 * T 
  //  - 0.0013268 * T * T + (T*T*T / 538841.0 ) - ( T*T*T*T / 65194000.0);

  // NEW VERSION
  /* mean length of moon */
  ls = 218.3164477 + 481267.88123421 * T 
    - 0.0015786 * T * T + (T*T*T / 538841.0 ) - ( T*T*T*T / 65194000.0);

  ls = fmod(ls , 360.0);
  if ( ls < 0 ) {
    ls += 360.0;
  }
  /* mean length of sun */
  l = 280.4665 + 36000.7698 * T ;
  l = fmod(l , 360.0);
  if ( l < 0 )
    l += 360.0;

  /* mean elongation of moo */
  //  d = 297.8502042 + 445267.1115168 * T 
  //  - 0.0016300 * T * T + ( T*T*T / 545868.0 ) - (T*T*T*T / 113065000.0);
  d = 297.8501921 + 445267.1114034 * T 
    - 0.0018849 * T * T + ( T*T*T / 545868.0 ) - (T*T*T*T / 113065000.0);

  d = fmod(d, 360.0);
  if ( d < 0 )
    d += 360.0;
  dmmf[0] = d;

  /* mean anomaly of sun */
  m = 357.5291092 + 35999.0502909 * T 
    - 0.0001536 * T * T + ( T*T*T / 24490000.0 );
  m = fmod(m, 360.0);
  if ( m < 0 ) 
    m += 360.0;
  dmmf[1] = m;

  /* mean anomaly of moon */
  //  ms = 134.9634114 + 477198.8676313 * T 
  //    + 0.0089970 * T * T + ( pow( T, 3) / 69699 ) - ( pow( T, 4) / 14712000);
  ms = 134.9633964 + 477198.8675055 * T 
    + 0.0087414 * T * T + ( T*T*T / 69699.0 ) - ( T*T*T*T / 14712000.0);
  ms = fmod(ms, 360.0);
  if ( ms < 0 )
    ms += 360.0;
  dmmf[2] = ms;

  /* mean distance of moon of rising knot */
  //  f = 93.2720993 + 483202.0175273 * T 
  //    - 0.0034029 * T * T - ( T*T*T / 3526000.0 ) - ( T*T*T*T / 863310000.0);
  f = 93.272095 + 483202.0175233 * T 
    - 0.0036539 * T * T - ( T*T*T / 3526000.0 ) - ( T*T*T*T / 863310000.0);
  f = fmod( f, 360.0);
  if ( f < 0 )
    f += 360.0;
  dmmf[3] = f;

  /* length of rising knot */
  Om = 125.04452 - 1934.136261 * T 
    + 0.0020708*T*T + T*T*T / 450000.0;
  Om = fmod( Om, 360.0);
  if ( Om < 0 )
    Om += 360.0;

  /* influence of venus, jupiter */
  a1 = 119.75 + 131.849 * T;
  a1 = fmod( a1, 360.0);
  
  a2 = 53.09  + 479264.290 * T;
  a2 = fmod (a2, 360.0);
  if ( a2 < 0 )
    a2 += 360.0;

  a3 = 313.45 + 481266.484 * T;
  a3 = fmod (a3, 360.0);
  if ( a3 < 0 )
    a3 += 360.0;
  
  /* excentricity of earth orbit around sun */
  e = 1.0 - ( 0.002516 * T ) - ( 0.0000074 * T * T );

  sbsum = 0.0;
  srsum = 0.0;
  slsum = 0.0;
  nterms = MOON_NCOEFF;

  for ( j = 0; j < nterms; j++) { 
    for ( i = 0; i < 4; i++) {
      tmpsb[i] = sb[j].arg[i];
      tmpsl[i] = slsr[j].arg[i];
    }
    dotsb   = dot( tmpsb, dmmf, 4);
    dotslsr = dot( tmpsl, dmmf, 4);
    sbres   = sin(deg2rad(dotsb)) * sb[j].sincoeff;
    slres   = sin(deg2rad(dotslsr)) * slsr[j].sincoeff;
    srres   = cos(deg2rad(dotslsr)) * slsr[j].coscoeff;

    if ( sb[j].emod == 0 ) 
      sbres = sbres * 1.0;
    else if ( sb[j].emod == 1 )
      sbres = sbres * e;
    else if ( sb[j].emod == 2)
      sbres = sbres * e * e;

    sbsum += sbres;
    srsum += srres;
    slsum += slres;
  }
  
  /* addtl. termes sum b */
  sbsum = sbsum - (2235.0 * sin(deg2rad(ls))) + ( 382.0 * sin(deg2rad(a3)));
  sbsum = sbsum
    + 175.0 * sin(deg2rad(a1 - f))
    + 175.0 * sin(deg2rad(a1 + f))
    + 127.0 * sin(deg2rad(ls - ms))
    - 115.0 * sin(deg2rad(ls + ms)) ;

  /* addtl terms sum l */
  slsum = slsum + 3958.0 * sin(deg2rad(a1))
    + 1962.0 * sin(deg2rad(ls -f ))
    + 318 * sin(deg2rad(a2));
  
  lambda = ls + ( slsum / 1000000.0 );
  beta   = sbsum / 1000000.0;
  delta  = 385000.56 + (srsum / 1000.0);
  sigma  = asin( 6378.14 / delta ) * 180.0 / ( PI ); 

  /* nutation */
  /* first approximation */
  dpsi = -17.20 * sin( deg2rad( Om)) 
    - 1.21 * sin( deg2rad( 2 * l ))
    - 0.23 * sin( deg2rad( 2 * ls)) 
    + 0.21 * sin( deg2rad( 2 * Om));
  dpsi = dpsi / 3600.0;

  lambda += dpsi;

  /* ecliptic calculations */
  deps = 9.20 * cos( deg2rad( Om))
    + 0.57 * cos( deg2rad( 2 * l))
    + 0.10 * cos( deg2rad( 2 * ls))
    - 0.09 * cos( deg2rad( 2 * Om));
  deps = deps /3600.0;

  U = T / 100.0;

  eps0 = - 4680.93 * U
    - 1.55 * U * U
    + 1999.25 * pow( U, 3)
    -   51.38 * pow( U, 4)
    -  249.67 * pow( U, 5)
    -   39.05 * pow( U, 6)
    +    7.12 * pow( U, 7)
    +   27.87 * pow( U, 8)
    +    5.79 * pow( U, 9)
    +    2.45 * pow( U, 10);
  eps0 = 23.439291 + eps0 /3600.0 ;
  eps = eps0 + deps;

  /* right ascension */
  Nom = sin( deg2rad( lambda)) * cos( deg2rad( eps)) 
    - tan( deg2rad( beta)) * sin( deg2rad( eps));
  Den = cos( deg2rad(lambda)); 
  ra = atan2( Nom, Den) * 180.0 / (PI) ;
  ra = fmod(ra, 360.0);
  while (ra < 0) ra += 360.0;
  while (ra >= 360) ra -= 360.0;

  /* declination */
  sindecl = sin( deg2rad( beta)) * cos( deg2rad( eps))
    + cos( deg2rad( beta)) * sin( deg2rad( eps)) * sin( deg2rad( lambda));
  decl = asin(sindecl) * 180.0 / (PI);

  MoonPar->ra = ra;
  MoonPar->decl = decl;

  return(0);
}

void moonpos(double jd, double *raout, double *decout) {

  OrbPar MoonCor;

  double t,ra,dec;

  t = ( jd - 2451545.0 ) / 36525.0;
  getmoonpos(t, &MoonCor);
  
  ra = (double) MoonCor.ra;
  dec = (double) MoonCor.decl;
  
  *raout=ra;
  *decout=dec;

}
