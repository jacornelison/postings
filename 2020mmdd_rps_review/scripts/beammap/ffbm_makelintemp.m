function T = ffbm_makelintemp(sig,ad);
% T = ffbm_makelintemp(sig,ad)
%
% Make a template equivalent to the deprojection templates, but in beam map
% space.  They are just derivatives of a circular gaussian.
%
% INPUT: sigma in degrees

% Make a Gaussian, since the T map
% Same x and y axes as standard maps

if ~exist('ad','var')
  x = -2:0.1:2;
  y = -2:0.1:2;
  y = -y;
  [X,Y] = meshgrid(x,y);
else
  [X,Y] = meshgrid(ad.t_val_deg{1},-ad.t_val_deg{2});
end

sigma_maj = sig;
sigma_min = sig;
alpha = 0;

beamparam = [1/(2*pi*sig^2) 0 0 sigma_maj sigma_min alpha];

% Diff relgain
T{1} = egauss2(beamparam,X,Y);

% Renormalize template
apeak = beamparam(1)/sum(sum(T{1}));
T{1}= T{1}./sum(sum(T{1}));

% Diff x pointing
T{2} = T{1}.*(X)/sig^2;

% Diff y pointing
T{3} = T{1}.*(Y)/sig^2;

% Differential beamwidth
T{4} = T{1}./sig^4.*(X.^2 + Y.^2 - 2*sig^2);

% Diff ellip: 
T{5} = (X.^2 - Y.^2).*T{1}./(sig^4);
T{6} = 2*X.*Y.*T{1}./sig^4;

return
