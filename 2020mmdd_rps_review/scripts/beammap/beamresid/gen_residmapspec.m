function [out beam split] = gen_residmapspec(residopt)
% [out beam split] = gen_residmapspec(residopt)
%
% Given a T map, differential Gaussian parameters OR actual beams, and a set
% of observing angles, determine the leaked Q/U maps and corresponding E/B
% spectra.
% 
% INPUTS: residopt struct containing fields
%           angles: vector of angles (radians) with which to observe the
%                   input T sky.  0 means the beam and map are simply
%                   convolved without rotation.  Can be as many as you
%                   want and the final map is the average of them.
%           type:   'image' or 'fourier' -- input T map and beam should
%                   have both in struct.  'image' is obviously much
%                   slower but seems to work better
%           ad_T:   ad struct corresponding to input T map
%           T:      temperature map with which to convolve the beam,
%                   should have fields 'im' and 'fp'
%           bintype: chooses binning in dumb_aps (e.g. 'bicep_norm')
%                    Or can send in n_bins to divide up ad equally
%
%        To generate an elliptical diff Gaussian with gen_egauss_par:
%           sigma:  nominal beamwidth, degrees
%           dparam: vector containing
%                   [dsig,    dx,     dy,     dp,       dc]
%                    degrees, arcmin, arcmin, unitless, unitless
%
%        To use actual beams, include the following
%           beam:   real beam
%                   should have fields 'im' and 'fp' for A, B, sum, diff
%           ad_beam: ad struct corresponding to beam, need not be the
%                    same as T (will interpolate later)
%           split:  optional, split map corresponding to noise estimate
%
% OUTPUTS:
%
%        out: struct containing Q/U maps and supfac-corrected E/B fourier
%             modes and binned EE/BB spectra
%        beam: beam used to convolve with T map (generated in this
%              function if just testing a Gaussian, or interpolated to T
%              map pixelization if we sent in a real beam)
%        split: same for split map if it exists

% Parse residopt
% These we have regardless of generating a beam or not
angles = residopt.angles;
type = residopt.type;
ad_T = residopt.ad_T; 
T = residopt.T;
bintype = residopt.bintype; 

% Load up or generate beam maps
if isfield(residopt,'sigma') 
  % Option 1: simple elliptical Gaussians - generate here
  sigma = residopt.sigma;
  dparam = residopt.dparam;
  beam = gen_egauss_pair(ad_T,sigma,dparam);
elseif isfield(residopt,'beam')
  % Option 2: real beam maps
  % Here we have 'beam' and 'ad_beam', get it into ad_T
  % Make sure to scale the difference!!!
  if isfield(residopt,'split')
    [beam split] = bm_interp(ad_T,residopt.ad_beam,...
                             residopt.beam,residopt.split,1);
  else
    beam = bm_interp(ad_T,residopt.ad_beam,residopt.beam,[],1);
  end
end

% Generate leakage maps
[Q U] = gen_QU_leak(ad_T,T,beam,angles,type);    
    
% reduc_makeaps works like so:
% ft = calc_map_fts(map,ad,weight)
% [l,Cs_l] = calcspec(ad,ft1,ft2,w,bintype)
% Above in principle I already calculated the map FTs for Q/U

% Turn Q/U into E/B
[u,v] = meshgrid(ad_T.u_val{1},ad_T.u_val{2});
[Ef_s Bf_s] = qu2eb(u,v,Q.fp,U.fp,[],[],'iau');

% Divide out beam suppression factor
% Beam sum power spectrum is included in 'beam'
Ef = Ef_s./beam.Sf;
Bf = Bf_s./beam.Sf;

[bc Csl_E] = aps_simple(ad_T,Ef,Ef,bintype);
[bc Csl_B] = aps_simple(ad_T,Bf,Bf,bintype);

out.l = bc;
out.EE = Csl_E;
out.BB = Csl_B;
out.Ef = Ef;
out.Bf = Bf;
out.Q = Q;
out.U = U;

% If split exists, do the same 
if isfield(residopt,'split')
  % Generate leakage maps
  [Q U] = gen_QU_leak(ad_T,T,split,angles,type);    
    
  % Turn Q/U into E/B
  [u,v] = meshgrid(ad_T.u_val{1},ad_T.u_val{2});
  [Ef_s Bf_s] = qu2eb(u,v,Q.fp,U.fp,[],[],'iau');

  % Use beam sum, not split sum for beam correction
  Ef = Ef_s./beam.Sf;
  Bf = Bf_s./beam.Sf;

  [bc Csl_E] = aps_simple(ad_T,Ef,Ef,bintype);
  [bc Csl_B] = aps_simple(ad_T,Bf,Bf,bintype);

  out.l_split = bc;
  out.EE_split = Csl_E;
  out.BB_split = Csl_B;
  out.Ef_split = Ef;
  out.Bf_split = Bf;
  out.Q_split = Q;
  out.U_split = U;
end
    
return

