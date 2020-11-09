function ffbm_param_manyfiles(experiment,year,rxNum)
% function ffbm_param_manyfiles(experiment,year,rxNum)
%
% Takes the output from beamfitter_par.m + output from 
% ffbm_cutting_pixels.m and summarizes the beam parameters 
% into a single file.
%
% The output from this file can be fed into
% ffbm_makebeamparamcsv.m to make beam parameter csv files.
%
% INPUTS
%         experiment:  'keck' or 'bicep2'
%         year:        2012, 2013, 2014
%         rxNum:       0/1/2/3/4 for Keck
%
% CLW 2014-06-16
% CLW 2014-06-18

% Get array info for specific year and load mask
switch experiment
  case 'bicep2'
    switch year
      case 2011
	[p ind] = get_array_info('20110201');
      case 2012
	[p ind] = get_array_info('20120201');
    end
    M = load(['bicep2_mask_' num2str(year) '/mask_b2.mat'])
  case 'keck'
    switch year
      case 2012
	[p ind] = get_array_info('20120201');
      case 2013
	[p ind] = get_array_info('20130201');
      case 2014
	[p ind] = get_array_info('20140201');
    end
    cutind = (p.rx == rxNum);
    p = structcut(p,cutind);
    ind = make_ind(p);
    p.theta2 = p.theta;
    p.theta = p.theta2 + p.drumangle;
    M = load(['mask_ffbm_' num2str(year) '/mask_rx' num2str(rxNum) '.mat'])
end

mask = M.mask;
maps = struct;
numFiles = length(mask)

for ii = 1:numFiles
  
  % Load up each beam map
  switch experiment
    case 'bicep2'
      bb = load(['beammaps/map_' char(mask(ii).time) '.mat']);
      aa = load(['beammaps/mapwin_' char(mask(ii).time) '.mat']);
      aa.dk = bb.dk;
      clear bb;
    case 'keck'
      aa = load(['beammaps/mapwin_' char(mask(ii).time) '_rx' num2str(rxNum) '.mat']);
  end

  % Calculate parameters with rotation
  [ss pp cc] = ffbm_findingparam(aa.A,p,aa.dk,1,experiment);
  sig(:,ii) = ss;
  e_p(:,ii) = pp;
  e_c(:,ii) = cc;

  % Calculate differential pointing
  clear g
  g = find_diffpoint_rot(aa.A,p,ind,aa.dk,experiment);
  x(:,ii) = g.x;
  y(:,ii) = g.y;

end

% Cut bad pixels
for ii = 1:numFiles
  % Make the master cut mask
  cutMask(:,ii) = mask(ii).poscut + ...
                  mask(ii).darksquid + ...
		  mask(ii).nofits + ...
		  mask(ii).badsig + ...
		  mask(ii).badellip + ...
		  mask(ii).nansigellip + ...
		  mask(ii).pkcut + ...
		  mask(ii).handcut;

  % NaN it out!	      
  for jj = 1:528
    if cutMask(jj,ii) ~= 0 % Bad measurement from flag
      sig(jj,ii) = NaN;
      e_p(jj,ii) = NaN;
      e_c(jj,ii) = NaN;
      x(jj,ii) = NaN;
      y(jj,ii) = NaN;
    end
  end

  %cut on sigma
  %there are some pixels with large spikes that i can't get rid of
  %takes out 1 measurement of 2 pixels in rx1
  %and 1 measurement of 9 pixels in rx2
  %and 1 measurement of 1 pixel in rx4
  %kk=find(sig(:,ii)>0.3 | sig(:,ii)<0.15)
  %sig(kk,ii)=NaN;
  %e_p(kk,ii)=NaN;
  %e_c(kk,ii)=NaN;
  %x(kk,ii)=NaN;
  %y(kk,ii)=NaN;
end

switch experiment
  case 'bicep2'
    save(['beamparam/beamparam_bicep2_' num2str(year)],'sig','e_p','e_c','x','y')  
  case 'keck'
    save(['beamparam_' num2str(year) '/beamparam_rx' num2str(rxNum)],'sig','e_p','e_c','x','y')
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = find_diffpoint_rot(A,p,ind,dk,experiment)

% Find rotation angle
switch experiment
  case 'bicep2'
    rotangle = -dk; % BICEP2's pointing reverses dk
  case 'keck'
    rotangle = p.drumangle(1) + dk;
end

rotangle = rotangle * pi/180; % Radians

% AB offset is A-B
d_az = A(2,ind.a) - A(2,ind.b);
d_el = A(3,ind.a) - A(3,ind.b);
d_el = -d_el; % The y axis points downwards, unfortunately

% So we rotate rotangle ccw
[theta,rho] = cart2pol(d_az,d_el);
thetanew = theta + rotangle;
[dx_p,dy_p] = pol2cart(thetanew,rho);

% Go from r and theta to x' and y'
[xcen,ycen] = pol2cart(p.theta*pi/180,p.r);

x(ind.a) = xcen(ind.a) + dx_p'/2;
y(ind.a) = ycen(ind.a) + dy_p'/2;

x(ind.b) = xcen(ind.b) - dx_p'/2;
y(ind.b) = ycen(ind.b) - dy_p'/2;

% [g.theta_ab,g.r_ab] = cart2pol(x,y);
% g.theta_ab = g.theta_ab*180/pi; % put things back into degrees
g.x = x;
g.y = y;

return