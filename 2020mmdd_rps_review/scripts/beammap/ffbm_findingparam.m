function [sig e_p e_c] = ffbm_findingparam(A,p,dk,rot,experiment)
% function [sig e_p e_c] = ffbm_findingparam(A,p,dk,rot,experiment)
%
% Takes the A matrix fit from  tes_analysis/util/normfit2d 
% (option to rotate them into the relevant x' and y' axis),
% and outputs sig, e_p, e_c in the x and y axis for which A is defined.
%
% x' and y' axes defined as tied to each focal plane.
% (The A matrix from  has the + y-axis downwards)
%
% Uses Stefan's function egauss_mmt2scp.m to convert ellipse params
%
% INPUTS
%         A:           Matrix from normfit2d, size(A) = [7,n_detectors]
%         p:           Only useful when we want to rotate into x' and y'
%                      From get_array_info (provides drumangle)
%         dk:          Dkangle for this map
%         rot:         0 (default, do not rotate the x/y axes to x'/y'
%                      1 to rotate to x'/y'
%         experiment:  'b2','keck','b3' (keck default)

if ~exist('experiment','var')
  experiment = 'keck';
end
if ~exist('rot','var')
  rot = 0;
  rotangle = 0;
end

if rot
  switch experiment
    case {'b2','bicep2'}
      rotangle = -dk + p.drumangle; % B2 pointing reverses apparent dk
    case {'keck','b3','bicep3'}
      rotangle = p.drumangle + dk;
  end
end

jj = 1;
num = length(A(1,:));
sig = NaN(1,num);
e_p = NaN(1,num);
e_c = NaN(1,num);
ellip = NaN(1,num);
angle = NaN(1,num);
angle2 = NaN(1,num);

for ii = 1:num
  if ~isnan(A(:,ii))
    
    varxy = A(6,ii) * A(5,ii) * A(4,ii);
    % 2x2 sigma matrix
    Sigma = [A(4,ii).^2 varxy ; varxy A(5,ii).^2];
    [V,D] = eigs(Sigma);
    
    % D is diag(sigma_major^2 sigma_minor^2)
    % y-axis is positive downwards, so theta is actually clockwise in matrix
    alpha = atan2(V(2,1),V(1,1)) * 180/pi;
    angle2(ii) = alpha;

    % Still need to account for the dk angle rotation:
    % rotangle = p.drumangle(1) + dk + 90; % CCW rotation in matrix space
    % This takes it to tile 1/2/3/4 the way we present it in postings.
    % We want x' to be approx along the line between tile 1 and 2
    %     and y' to be approx along the line between tile 1 and tile 4...
    % So rotangle = p.drumangle(1) + dk;
    % dk angle is rotated CW when viewed from inside the sphere onto the sky
    % dk angle is a clockwise rotation
    % So to remove dk angle, you have to rotate CCW.
    % -> we rotate dk CCW, and drumangle CCW.
    % alpha is clockwise in the matrix
    if rot
      alpha = -alpha; % Make alpha CCW
      %alpha = alpha + rotangle; % Rotate CCW to ccount for dk + drumangle angle
      alpha = alpha + rotangle(ii);
      angle(ii) = alpha;
    else
      angle(ii) = alpha;
    end
    
    fwhm_maj = 2 * sqrt(2*log(2)) * sqrt(D(1,1));
    fwhm_min = 2 * sqrt(2*log(2)) * sqrt(D(2,2));
    [sig(ii),e_c(ii),e_p(ii)] = egauss2_mmt2scp(fwhm_maj,fwhm_min,angle(ii));
    ellip(ii) = (sqrt(D(1,1)) - sqrt(D(2,2))) / sig(ii)/2;
    
    if ~isreal(sig(ii))
      % Get rid of imaginary values
      sig(ii) = NaN; 
      disp(sprintf('Bad Sigma det: %d',ii))
      badsig(jj) = ii;
      jj = jj + 1;
    end
    if ~isreal(ellip(ii))
      ellip(ii) = NaN;
      disp(sprintf('Bad Ellipticity det: %d ',ii))
      badsig(jj) = ii;
      jj = jj + 1;
    end
    
  end
end


return