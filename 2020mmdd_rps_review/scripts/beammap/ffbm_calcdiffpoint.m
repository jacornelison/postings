function [x y] = ffbm_calcdiffpoint(A,p,ind,dk,expt,rot)
% [x y] = ffbm_calcdiffpoint(A,p,ind,dk,experiment,rot)
%
% Takes the A matrix fit from tes_analysis/util/normfit2d, rotates into x'
% and y', and outputs the beam positions.
% Function name is a bit misleading since you need to difference x and y to
% get differential pointing.
% Pulled from ffbm_param_manyfiles
%
% INPUTS
%         A:           Matrix from normfit2d, size(A) = [7,n_detectors]
%         p:           From get_array_info (provides drumangle)
%         ind:         From get_array_info (provides A/B indices)
%         dk:          Dkangle for this map
%         expt:        'b2','keck','b3' (keck default)
%         rot:         = 1 (default) rotate to x'/y'
%                      = 0 no rotation applied

if ~exist('rot','var')
  rot = 1;
end

% Find rotation angle
if rot
  switch expt
    case {'b2'}
      rotangle = -dk + p.drumangle ; % BICEP2's pointing reverses dk
    case {'keck','b3','bicep3'}
      rotangle = p.drumangle + dk;
  end
else
  rotangle = zeros(size(p.drumangle));
end

rotangle = rotangle * pi/180; % Radians
rotangle = rotangle(ind.a); % Get the right dimensions (per-pair)

% AB offset is A-B
d_az = A(2,ind.a) - A(2,ind.b);
d_el = A(3,ind.a) - A(3,ind.b);
d_el = -d_el; % The y axis points downwards, unfortunately
              % This is effectively undoing the mirror inversion
              % See 20150714_bmaxes posting

% So we rotate rotangle ccw
[theta,rho] = cart2pol(d_az,d_el);
thetanew = theta + rotangle';
[dx_p,dy_p] = pol2cart(thetanew,rho);

% Go from r and theta to x' and y'
[xcen,ycen] = pol2cart(p.theta*pi/180,p.r);

x = NaN(length(ind.e),1);
y = NaN(length(ind.e),1);

x(ind.a) = xcen(ind.a) + dx_p'/2;
y(ind.a) = ycen(ind.a) + dy_p'/2;

x(ind.b) = xcen(ind.b) - dx_p'/2;
y(ind.b) = ycen(ind.b) - dy_p'/2;

% [g.theta_ab,g.r_ab] = cart2pol(x,y);
% g.theta_ab = g.theta_ab*180/pi; % put things back into degrees
%g.x = x;
%g.y = y;

return