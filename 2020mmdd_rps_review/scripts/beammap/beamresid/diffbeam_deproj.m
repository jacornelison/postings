function dpbeam = diffbeam_deproj(beam,ad,sigma,dp)
% dpbeam = diffbeam_deproj(beam,ad,sigma,dp)
%
% Take the beam and deproject it in map space
% Only works on image-space difference beam
% Only works with 4-character dp label for now
%
% INPUTS
%         beam: real beam
%               Should have fields 'i' and 'f' for A/B/sum/diff
%               (but this only works on Di)
%         ad:   standard struct for beam map
%         sigma: nominal beamwidth (degrees) for generating 
%                GAUSSIAN deprojection template
%         dp:    string corresponding to modes that should be
%                deprojected, e.g. 1100 for relgain/diffpoint
%                Does not subtract yet
% 
% OUTPUTS
%         dpbeam: beam struct with beam.Di deprojected

% Just like deprojection in the pipeline, we have to simultaneously
% deproject only the modes we want (not all at once)
modes = [];
if strcmp(dp(1),'1')
  modes = [modes 1];
end
if strcmp(dp(2),'1')
  modes = [modes 2 3];
end
if strcmp(dp(3),'1')
  modes = [modes 4];
end
if strcmp(dp(4),'1')
  modes = [modes 5 6];
end

[coeff T] = ffbm_deprojmap(beam.Di,ad,sigma,modes);

toremove = 0;

% Relgain
if strcmp(dp(1),'1')
  toremove = toremove + coeff(1).*T{1};
end

% Diff point
if strcmp(dp(2),'1')
  toremove = toremove + coeff(2).*T{2} + coeff(3).*T{3};
end

% Diff beamwidth
if strcmp(dp(3),'1')
  toremove = toremove + coeff(4).*T{4};
end

% Diff ellip
if strcmp(dp(4),'1')
  toremove = toremove + coeff(5).*T{5} + coeff(6).*T{6};
end

dpbeam = beam;
dpbeam.Di = dpbeam.Di - toremove;
  
return