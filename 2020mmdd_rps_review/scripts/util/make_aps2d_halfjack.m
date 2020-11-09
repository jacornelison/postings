function [ad,aps2d,map,nonjackmap]=make_aps2d_halfjack(m,map,nonjackmap,pure_b, pad, realimag)
% [ad,aps2d,map]=make_aps2d_halfjack(m,map,pure_b)

% stolen from make_aps2d
% takes 2 maps as the input (a jackknifed map, and one of the 2 original maps)
% returns the "half jack" cross spec
% a normal jackknife cross spec will return (T1-T2)*conj(E1-E2)
% while a halfjack returns T1*conj(E1-E2)

if ~exist('pad','var')
  pad=1;
end

if ~exist('realimag','var')
  realimag=0;
end


if ~exist('pure_b','var')
  pure_b='normal';
end


% needs to be calibrated, var maps need to be smoothed.
% pad the maps
if(pad)
  disp('pad')
  tic
  pp=2^(nextpow2(max([m.nx,m.ny]))+1);
  [m,map]=pad_map(m,map,pp);
  [m,nonjackmap]=pad_map(m,nonjackmap,pp);
  toc
end

% calc axis data
ad=calc_ad2([m.xdos,m.ydos],[m.nx,m.ny]);

disp('fft')
tic
% calculate the fourier components of maps
map=calc_map_fouriercomp(map,ad,pure_b);
nonjackmap=calc_map_fouriercomp(nonjackmap,ad,pure_b);
toc

% calc 100/150GHz 2d aps
sf=prod(ad.del_u);
disp('calc 100/150GHz 2d aps');
for j=1:numel(map)
  if isfield(map(j),'Q')    
    aps2d(j).TE=specplot(nonjackmap(j).Tft, map(j).Eft,sf,realimag);
    aps2d(j).TB=specplot(nonjackmap(j).Tft, map(j).Bft,sf,realimag);
    aps2d(j).EB=specplot(nonjackmap(j).Bft, map(j).Bft,sf,realimag);
    aps2d(j).TQ=specplot(nonjackmap(j).Tft, map(j).Qft,sf,realimag);
    aps2d(j).TU=specplot(nonjackmap(j).Tft, map(j).Uft,sf,realimag);
  end
end

% mult to flat in l(l+1)C_l/2pi
ad.l_r=2*pi*ad.u_r;
sf=ad.l_r.*(ad.l_r+1)/(2*pi);


 %%%%%%%%%%%%%%%%%%%%
function spec_out=specplot(spec_in_A,spec_in_B,sf,realimag)
if(~realimag)
spec_out=sf*real(spec_in_A.*conj(spec_in_B));
else  
width=size(spec_in_A,2);
imagspec=sf*imag(spec_in_A).*imag(spec_in_B);
realspec=sf*real(spec_in_A).*real(spec_in_B);
spec_out(:,1:floor(width/2))=imagspec(:,1:floor(width/2));
spec_out(:,floor(width/2)+1:width)=realspec(:,floor(width/2)+1:end);
end

return
