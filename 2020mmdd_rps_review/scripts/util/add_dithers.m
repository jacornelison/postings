function p=add_dithers(p,pp,rlz,simoptrlz)
% p=add_dithers(p,pp,rlz,simopt.rlz)
%
% Modify p structure with dithers in pp structure for the correct realization.

% If doing rlz=0, use first dither. This is a bit of a kludge to allow dithering of a
% parameter when using a WMAP/Planck map as an input.
if rlz==0
  rlz=1;
  simoptrlz=1;
end

% Are the dither parameters Nchan x 1 or Nchan x Nrlz? If the former, always pick out
% the first column. If the latter, pick out the correct column as specified in rlz.
fld=fieldnames(pp);
rlzind=[];
for i=1:numel(fld)
  x=getfield(pp,fld{i});
  if size(x,2)==numel(simoptrlz)
    % this is a dithered parameter
    rlzind=setfield(rlzind,fld{i},find(simoptrlz==rlz));
  elseif size(x,2)==1
    rlzind=setfield(rlzind,fld{i},1);
  else
    % Something is wrong
    error('Dither parameter size must be size Nchan x max(simopt.rlz) or Nchan x 1');
  end
end

% dither epsilon 
p.epsilon=p.epsilon+pp.epsilon(:,rlzind.epsilon);
% dither chi
p.chi=p.chi+pp.chi(:,rlzind.chi);

% dither beam centers
[ra_off_dos,dec_off]=pol2cart(p.theta*pi/180,p.r);
ra_off_dos=ra_off_dos+pp.ra_off_dos(:,rlzind.ra_off_dos);
dec_off=dec_off+pp.dec_off(:,rlzind.dec_off);
[p.theta,p.r]=cart2pol(ra_off_dos,dec_off);
p.theta=p.theta*180/pi;

%  % dither beam widths % prepare for removal (STF)
%  p.fwhm_maj=p.fwhm_maj+pp.fwhm_maj(:,rlzind.fwhm_maj);
%  p.fwhm_min=p.fwhm_min+pp.fwhm_min(:,rlzind.fwhm_min);
%  % dither beam angle
%  p.alpha=p.alpha+pp.alpha(:,rlzind.alpha);

% dither ellipticity and beam width using sigma,c and p:
% first fetch what is the status
[sigma,elc,elp] =  egauss2_mmt2scp(p.fwhm_maj,p.fwhm_min,p.alpha+p.theta);

% then dither ellipticity 
elp=(elp+pp.p(:,rlzind.p));
elc=(elc+pp.c(:,rlzind.c));

% now dither beamwidths: that is special as it has both an additive and a multiplicatice term (fist case)
% However, if dithers are coming from individual uncertainty on a channel by channel
% basis, the multiplicative term does not apply (second case) 
if (size(pp.sigma,3)==2)
  sigma=(sigma+pp.sigma(:,rlzind.sigma,1)).*pp.sigma(:,rlzind.sigma,2);
else
  sigma=(sigma+pp.sigma(:,rlzind.sigma));
end
% and convert back into fwhm_maj, fwhm_min, alpha
[p.fwhm_maj,p.fwhm_min,p.alpha] =  egauss2_scp2mmt(sigma,elc,elp);
p.alpha = p.alpha - p.theta;

% copy in random abgain
p.abgain=pp.abgain(:,rlzind.abgain);


return

