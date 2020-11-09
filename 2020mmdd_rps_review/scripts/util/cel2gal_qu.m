%Converted from Cynthia's C code "cel2gal_qu.c" in bicep_point.c
%JET 20111011
%
%Input: (qi, ui) = input maps
%       (ra, dec)= coordinate of input maps in cel

%Ouput: (qo, uo) = output maps
%       (l,b)    = output coord in gal

function [qo,uo, l, b]=cel2gal_qu(ra,dec, qi,ui)

    %/* J2000 */
    ra_gal_pole = (12.+51.4/60.)*pi/12.;
    dec_gal_pole = (27. + 7./60.)*pi/180.;
    ra_gal_origin = (17. + 45.6/60.)*pi/12.;
    dec_gal_origin=-(28.0+ 56./60.)*pi/180.;
    
    %/* Preliminary calculations for coordinate transformation */
    %use euler instead of coord()
      [l_cel_pole, b_cel_pole]=euler(0, 90, 1);
      l_cel_pole=l_cel_pole.*pi/180;
      b_cel_pole=b_cel_pole.*pi/180;
      [l_cel_origin, b_cel_origin]=euler(0,0,1);
      l_cel_origin=l_cel_origin.*pi/180;
      b_cel_origin=b_cel_origin.*pi/180;

    %use euler instead of cel2gal to find coords of new map
    [l,b]=euler(ra, dec,1);

    %/* beta = intersection angle between lines of constant RA and
    %   constant l running through point of interest */
    qtmp = qi; utmp = ui;
    p = sqrt(qtmp.*qtmp+utmp.*utmp);
    a = 0.5*atan2(utmp,qtmp);

    aa = pi/2-dec.*pi/180;
    cc = pi/2-b.*pi/180;
    
    %Law of Cosines
    betan=cos(pi/2-b_cel_pole) - cos(cc).*cos(aa);
    betad=sin(cc).*sin(aa);
    beta=betan./betad;
    beta = acos(beta);

 %/* Sign of beta */
   ff=ra*pi/180-ra_gal_pole;
   f=find(ff<0);
     	beta(f) = -1.0*beta(f);
 
   %now do subtraction
    a = a-beta;

    %only real part
    a=real(a);

    %/* New Q and U */
    qo = p.*cos(2*a);
    uo = p.*sin(2*a);
   
 return