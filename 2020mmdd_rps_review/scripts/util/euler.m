function [ao, bo] = euler(ai,bi,select, FK4)
% [ao, bo] = euler(ai,bi,select, FK4)
%
% in/out in DEGREES
%
% e.g. [l,b]=euler(ra,dec,1,0)
%      [ra,dec]=euler(l,b,2,0)

% NAME:
%     EULER
% NOTE: HACKED by JMK 1/13/00 to translate to MATLAB from IDL
%
% PURPOSE:
%     Transform between Galactic, celestial, and ecliptic coordinates.
% EXPLANATION:
%     Use the procedure ASTRO to use this routine interactively
% 
% CALLING SEQUENCE:
%      EULER, AI, BI, AO, BO, [ SELECT, /FK4 ] 
%
% INPUTS:
%       AI - Input Longitude in DEGREES, scalar or vector.  If only two 
%               parameters are supplied, then  AI and BI will be modified to 
%               contain the output longitude and latitude.
%       BI - Input Latitude in DEGREES
%
% OPTIONAL INPUT:
%       SELECT - Integer (1-6) specifying type of coordinate transformation.  
%
%      SELECT   From          To        |   SELECT      From            To
%       1     RA-Dec (2000)  Galactic   |     4       Ecliptic      RA-Dec    
%       2     Galactic       RA-DEC     |     5       Ecliptic      Galactic  
%       3     RA-Dec         Ecliptic   |     6       Galactic      Ecliptic  
%
%      If omitted, program will prompt for the value of SELECT
%      Celestial coordinates (RA, Dec) should be given in equinox J2000 
%      unless the /FK4 keyword is set.
% OUTPUTS:
%       AO - Output Longitude in DEGREES
%       BO - Output Latitude in DEGREES
%
% INPUT KEYWORD:
%       /FK4 - If this keyword is set and non-zero, then input and output 
%             celestial and ecliptic coordinates should be given in equinox 
%             B1950.
%
% NOTES:
%       EULER was changed in December 1998 to use J2000 coordinates as the 
%       default, ** and may be incompatible with earlier versions***.
% REVISION HISTORY:
%       Written W. Landsman,  February 1987
%       Adapted from Fortran by Daryl Yentis NRL
%       Converted to IDL V5.0   W. Landsman   September 1997
%       Made J2000 the default, added /FK4 keyword  W. Landsman December 1998
%-
% On_error,2

% npar = N_params()
if (nargin < 2)
   disp('Syntax -  [A0, B0] = euler(AI,BI,SELECT, FK4)')
   disp('    AI,BI - Input longitude,latitude in degrees')
   disp('    AO,BO - Output longitude, latitude in degrees')
   disp('    SELECT - Scalar (1-6) specifying transformation type (default = 1)')
   disp('    FK4 - Boolean indicating equinox B1950 instead of J2000 (default = 0)')
   return
end

twopi   =   2 * pi;
fourpi  =   4 * pi;
deg_to_rad = 180 / pi;

%   J2000 coordinate conversions are based on the following constants
%  eps = 23.4392911111d              Obliquity of the ecliptic
%  alphaG = 192.85948d               Right Ascension of Galactic North Pole
%  deltaG = 27.12825d                Declination of Galactic North Pole
%  lomega = 32.93192d                Galactic longitude of celestial equator  
%  alphaE = 180.02322d              Ecliptic longitude of Galactic North Pole
%  deltaE = 29.811438523d            Ecliptic latitude of Galactic North Pole
%  Eomega  = 6.3839743d              Galactic longitude of ecliptic equator              

if exist('select','var')~=1
    select = 1;
end

if exist('FK4','var')~=1
    FK4 = 0;
end



if FK4
  equinox = '(B1950)';
  psi   = [ 0.57595865315, 4.9261918136, ...
            0.00000000000, 0.0000000000, ... 
            0.11129056012, 4.7005372834];
  stheta =[ 0.88781538514,-0.88781538514, ...
            0.39788119938,-0.39788119938, ...
            0.86766174755,-0.86766174755];   
  ctheta =[ 0.46019978478, 0.46019978478, ...
            0.91743694670, 0.91743694670, ...
            0.49715499774, 0.49715499774];
   phi  = [ 4.9261918136,  0.57595865315, ...
            0.0000000000, 0.00000000000, ...
	    4.7005372834, 0.11129056012];

else 

  equinox = '(J2000)';
  psi   = [ 0.57477043300, 4.9368292465,  ...
            0.00000000000, 0.0000000000,  ... 
            0.11142137093, 4.71279419371];    
  stheta =[ 0.88998808748,-0.88998808748, ...
            0.39777715593,-0.39777715593, ...
            0.86766622025,-0.86766622025];   
  ctheta =[ 0.45598377618, 0.45598377618, ...
            0.91748206207, 0.91748206207, ...
            0.49714719172, 0.49714719172];   
   phi  = [ 4.9368292465,  0.57477043300, ...
            0.0000000000, 0.00000000000, ...
            4.71279419371, 0.11142137093];

end


% if npar LT 5 then begin
%        print,' '
%        print,' 1 RA-DEC ' + equinox + ' to Galactic
%        print,' 2 Galactic       to RA-DEC' + equinox
%        print,' 3 RA-DEC ' + equinox + ' to Ecliptic'
%        print,' 4 Ecliptic       to RA-DEC' + equinox
%        print,' 5 Ecliptic       to Galactic'
%        print,' 6 Galactic       to Ecliptic'
%        select = 0
%        read,'Enter selection: ',select
% endif

 i  = select;                        % no matlab offset ; IDL offset
 a  = ai./deg_to_rad - phi(i);
 b = bi./deg_to_rad;
 sb = sin(b);
 cb = cos(b);
 cbsa = cb .* sin(a);
 b  = -stheta(i) .* cbsa + ctheta(i) .* sb;
 nanb = isnan(b);
 b  = min(b,1);
 b(nanb) = NaN;
 bo    = asin(b) .* deg_to_rad;

 a =  atan2( ctheta(i) .* cbsa + stheta(i) .* sb,  cb .* cos(a) );
 ao = (mod((a+psi(i)+fourpi), twopi)) * deg_to_rad;


% if ( npar EQ 2 ) then begin
%	ai = ao & bi=bo
% endif

 return
