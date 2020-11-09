function [be,n]=get_bins(bintype)
%function [be,n]-get_bins(bintype)
%function to return the l bin edges given a bintype
%should be called instead of defining band edges 

 switch bintype
   
  case 'bicep_fine'
   r=4;
   n=35/r; be=20:n:581;  be=[linspace(0,be(1),r+1),be(2:end)];
   

  case 'bicep_norm'
   r=1;
   n=35/r; be=20:n:581;  be=[linspace(0,be(1),r+1),be(2:end)];
   
  case 'bpwf'
   % same as quad fine but extends to 2750 to catch trailing edge of
   % last bin
   r=4;
   n=81/r; be=1.5:n:2751;
   
  case 'kappa'
   r=1;
   n=35/r; be=20:n:701;  be=[linspace(0,be(1),r+1),be(2:end)];
 
 case 'phi' 
   %bins in phi rather then ell
   n=(pi/2)*0.1;
   be=(pi/2)*(-1.05:0.1:1.05);
   
  otherwise
   error('unknown bintype')
   
 end
 
 
 return
 
