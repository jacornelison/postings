function ffbm_rotatemaps(compopt)
% ffbm_rotatemaps(times,experiment,rxNum,dk,choplet,mask)
% 
% Rotate FFBMs to the right orientation in order to make composite maps.
% We have two potential orientations:
% reduc_makesim wants maps defined such that the beams for the entire
% array projected onto the sky at dk = 0, i.e. with the drum angle built
% into the maps.
% However, for beam parameter comparison we also want to look at maps in
% x'/y' 
% PLEASE READ THE CODE BEFORE ATTEMPTING TO PLOT
%
% Centering on the A/B centroid or per-detector and masking should be done
% before this step.
%
% INPUTS (should all be sent in with compopt)
%
%  expt:          'bicep2','keck','bicep3'
%  year:          bicep2: 2012
%                 keck: 2012-2016
%                 bicep3: 2015
%
% OPTIONAL INPUTS
%
%  mapsize:       map size in deg (2 default), used for map filename
%  maskeddir:     directory from which to load masked maps
%                 (default beammaps/maps_masked)
%                 if needed can point somewhere else (if don't want masks)
%  maskedfile:    file in maskeddir to load
%                 (default ffbm_year_allcomp_masked)
%  rotateddir:    directory in which to save rotated maps
%                 (default beammaps/maps_rotated)
%  rotatedfile:   file in rotateddir to save
%                 (default ffbm_year_allcomp_rotated)
%  coord:         Coordinate system to rotate to:
%                 'dk0' (default) is the standard input to reduc_makesim
%                 'xpyp' for normal visualization (x'/y')
%  suffix:        optional string to add to map name, '' default

% Parse compopt
expt = compopt.expt;
year = compopt.year;
if ~isfield(compopt,'mapsize');
  sizestr = '2deg';
else
  sizestr = [strrep(num2str(compopt.mapsize),'.','p') 'deg'];
end
if ~isfield(compopt,'suffix')
  suffix = '';
else
  suffix = compopt.suffix;
end
if ~isfield(compopt,'maskeddir')
  compopt.maskeddir = 'beammaps/maps_masked';
  maskeddir = compopt.maskeddir;
else
  maskeddir = compopt.maskeddir;
end
if ~isfield(compopt,'maskedfile')
  compopt.maskedfile = ['ffbm_' num2str(year) '_allcomp_masked'...
	'_' sizestr suffix];
  maskedfile = compopt.maskedfile;
else
  maskedfile = compopt.maskedfile;
end
if ~isfield(compopt,'rotateddir')
  compopt.rotateddir = 'beammaps/maps_rotated';
  rotateddir = compopt.rotateddir;
else
  rotateddir = compopt.rotateddir;
end
if ~isfield(compopt,'rotatedfile')
  compopt.rotatedfile = ['ffbm_' num2str(year) '_allcomp_rotated'...
	'_' sizestr suffix];
  rotatedfile = compopt.rotatedfile;
else
  rotatedfile = compopt.rotatedfile;
end
if ~isfield(compopt,'coord')
  compopt.coord = 'dk0';
  coord = compopt.coord;
else
  coord = compopt.coord;
end

[p ind] = get_array_info([num2str(year) '0201']);

% Load up map file
filename = [maskeddir '/' maskedfile];
load(filename);

% Go through all light pairs and do yo thang
for ii = 1:length(ind.la)
  if mod(ii,10) == 0
    disp(['Rotating pair ' num2str(ii) ' / ' num2str(length(ind.la))]);
  end
  for jj = 1:length(comp.bm.number)

    % Component maps here have dimension 1 = increasing ap el
    %                          dimension 2 = increasing ap az
    % In ffbm_plotCompMaps this is rendered as
    % imagesc(comp.az_ap,comp.el_ap,map); set(gca,'ydir','normal') or
    % imagesc(comp.ad.t_val_deg{1},comp.ad.t_val_deg{2},map);
    %   set(gca,'ydir','normal)
    % where ad.t_val_deg{1} is increasing apparent az - mean(az)
    %       ad.t_val_deg{2} is increasing apparent el - mean(el)

    map_a = comp.map{ind.la(ii)}.component{jj};
    map_b = comp.map{ind.lb(ii)}.component{jj};
    
    % Rotate!  
    
    switch coord
      case 'dk0'
        % For full discussion, see again the 20150714_bmaxes posting
        % The sim currently (2016-05-02) wants the beam map to be in
        % instantaneous projection onto the sky of the beam where 
        % dimension 1 = decreasing elevation and dimension 2 = increasing
        % azimuth.  If dk = drumangle = 0, the beam map should be such that
        % dimension 1 = decreasing y' and dimension 2 = increasing x',
        % and if drumangle = 90, 1 = increasing x' and 2 = increasing
        % y'.  With the drum angle baked in, this means we simply rotate
        % by the dk angle of the beam map.
        rotangle = comp.bm.dk(jj); 
        % Fig. 10 shows that this places dimension 1 = decreasing
        % elevation and dimension 2 = increasing azimuth
      case 'xpyp'
        % All right bro
        % Due to the mirror flip, we need to reverse the parity.
        map_a = flipud(map_a);
        map_b = flipud(map_b);
        % Now rotate by the right angle:
        % We want dimension 1 = increasing x prime
        %         dimension 2 = decreasing y prime
        switch expt
          case 'bicep2'
            % CHECK ME
            rotangle = comp.bm.dk(jj) + p.drumangle(ind.la(ii)); 
          case 'keck'
            % We find empirically that Keck needs this extra 90 degrees.
            % Oh well.
            %rotangle = comp.bm.dk(jj) + p.drumangle(ind.la(ii)) + 90; 
            rotangle = -(90 + comp.bm.dk(jj) + p.drumangle(ind.la(ii)));
	  case 'bicep3'
	    % CHECK ME
	    rotangle = -90 + (comp.bm.dk(jj) + p.drumangle(ind.la(ii)));
        end
    end
    
    comp.map{ind.la(ii)}.component{jj} = ...
	imrotate(map_a,rotangle,'bilinear','crop');
    comp.map{ind.lb(ii)}.component{jj} = ...
	imrotate(map_b,rotangle,'bilinear','crop');
    
  end % beam map loop
end % Light A loop

% We now define what the ad axis means for our two coordinate systems

switch coord
  case 'dk0'
    % As explained above, the map array has
    % dimension 1 = decreasing elevation (increasing RA at dk=0)
    % dimension 2 = increasing azimuth   (increasing dec at dk=0)
    % So here we define 
    % comp.ad.t_val_deg{1} = decreasing elevation,
    % comp.ad.t_val_deg{2} = increasing azimuth 
    % but since the raw numbers of this array are actually increasing
    % as-is, we flip it here:  
    comp.ad.t_val_deg{1} = fliplr(comp.ad.t_val_deg{1});
    % BE VERY CAREFUL WHEN PLOTTING THIS WITH IMAGESC since it will
    % re-orient the array if you're not careful
    %   imagesc(ad.t_val_deg{2},ad.t_val_deg{1},map);
    %   set(gca,'ydir','normal'); xlabel('az'); ylabel('el')
    % should get you the right orientation if you want to plot this
  case 'xpyp'
    % For x prime y prime, we see that the above map array has
    % dimension 1 = increasing x prime
    % dimension 2 = decreasing y prime
    % We now define 
    % comp.ad.t_val_deg{1} = increasing x prime 
    % comp.ad.t_val_deg{2} = decreasing y prime
    % but again, since the raw numbers of the array are actually
    % increasing as-is, we flip it here:
    comp.ad.t_val_deg{2} = fliplr(comp.ad.t_val_deg{2});
    % So to plot our standard-parity beam plot, we do:
    %   imagesc(ad.t_val_deg{2},ad.t_val_deg{1},map);
    %   set(gca,'ydir','normal'); set(gca,'xdir','reverse')
    %   xlabel('y prime'); ylabel('x prime')
    % to get the numbers running backwards along the bottom
end

comp.coord = coord; % Save for plotting later
comp.compopt = compopt;

savename = [rotateddir '/' rotatedfile '_' coord];
save(savename,'comp','-v7.3');

disp(['Saved rotated component file: ' savename]);

return
%{ 

% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      experiment:        'bicep2','keck'
%      mapsourcetype:     '','uberchopper'
%      year:              2012,2013,2014 (keck only)
%      rxNum:             keck only
%      choplet:           0 (default), normal maps
%                         1, look for map_cl_
%      mask:              1 (default), mask ground/SPT
%                         0, no masking
%
% EXAMPLES
%
%      ffbm_rotatemaps('keck','',2012,1,-122,1)
%      ffbm_rotatemaps('keck','',2012,1,-50,1)
%      ffbm_rotatemaps('keck','',2012,1,22,1)
%      ffbm_rotatemaps('keck','',2012,1,-86,1)
%      ffbm_rotatemaps('keck','',2012,1,-14,1)
%
% CLW 2014-05-25
% KSK 2014-08-11
%

%[time dkangle sched] = ffbm_findbmruntimes(experiment,mapsourcetype,year,rxNum);  

% What times do we want?  Match both dk and schedule
%whichtime = find(dkangle == dk & sched == sc);

for ii = 1:length(times)
  clear map
  switch experiment
    case 'bicep2'
      f = load(['beammaps/maps/map_' char(times{ii}) '.mat']);
      w = load(['beammaps/maps/mapwin_' char(times{ii}) '.mat']);
    case 'keck'
      if choplet
	f = load(['beammaps/maps/map_cl_' char(times{ii}) '_rx' ...
	      num2str(rxNum) '.mat']);
	w = load(['beammaps/maps/mapwin_cl_' char(times{ii}) '_rx' ...
	      num2str(rxNum) '.mat']);
      else
	f = load(['beammaps/maps/map_' char(times{ii}) '_rx' ...
	      num2str(rxNum) '.mat']);
	w = load(['beammaps/maps/mapwin_' char(times{ii}) '_rx' ...
	      num2str(rxNum) '.mat']);
      end
  end
  
  % Find the median fit coordinates
  elcenter = nanmedian(w.A(3,:));  
  azcenter = nanmedian(w.A(2,:));
  
  f.mapmasked = f.map;
  % Mask out the ground completely
  if mask
    xlimcalc = f.x_bin;
    ylimcalc = f.y_bin < elcenter - 1.5;
    [calcmaskaz calcmaskel] = meshgrid(xlimcalc,ylimcalc);
    groundmask = calcmaskaz & calcmaskel;
    groundmask = repmat(groundmask,[1,1,528]);
    %f.mapmasked = f.map;
    f.mapmasked(groundmask) = NaN;
  end
    
  % Find dk rotation and if keck, mask out SPT
  switch experiment
    case 'bicep2'
      dkrot = -double(f.dk); % B2's pointing reverses dk angle      
    case 'keck'
      if mask
	xlimcalc = f.x_bin < azcenter - 2.5 & f.x_bin > azcenter - 7.5;
	ylimcalc = f.y_bin < elcenter + 1;
	[calcmaskaz calcmaskel] = meshgrid(xlimcalc,ylimcalc);
	sptmask = calcmaskaz & calcmaskel;
	sptmask = repmat(sptmask,[1,1,528]);
	%f.mapmasked = f.mapmasked;
	f.mapmasked(sptmask) = NaN;
      end
      dkrot = double(f.dk);
  end

  tic
  for jj = 1:528
    A = f.mapmasked(:,:,jj);
    map(:,:,jj) = imrotate(f.mapmasked(:,:,jj),dkrot,'bicubic');
  end
  toc
  ww{ii}.map = map;
  
  % This map is projected onto the sky at dk = 0, looking from the inside of
  % the sphere (when using imagesc(map)).  Compromise with y-axis flip:
  
  x_res = f.x_bin(2) - f.x_bin(1);
  y_res = f.y_bin(2) - f.y_bin(1);
  [ysize xsize] = size(ww{ii}.map(:,:,1));
  ww{ii}.x_bin = ([1:xsize] - xsize/2)*x_res;
  y_bin = ([1:ysize] - ysize/2)*y_res;
  ww{ii}.y_bin = fliplr(y_bin);
  
  switch experiment
    case 'bicep2'
      cent.x = 0;
      cent.y = 0;
    case 'keck'
      % azcenter and elcenter gives the median location of the beam center
      % before rotation.  We rotate CCW by dkrot about the center of the map.
      % The center is:
      xcentb = (f.x_bin(end) + f.x_bin(1))/2;
      ycentb = (f.y_bin(end) + f.y_bin(1))/2;
      % REMEMBER: y-axis increases downwards!
      % -> Theta here is CW, so to rotate CCW, we subtract
      [THETA,RHO] = cart2pol(azcenter - xcentb,elcenter - ycentb);
      [X,Y] = pol2cart(THETA - dkrot*pi/180,RHO);
      cent.x = X;
      % But we've flipped the y axis for the rotated map.  This really only
      % works because the center is 0.
      cent.y = -Y;
  end
    
  % Find beam centers
  for jj = 1:528
    
    jj
    if(1)
      time = times{ii}
      fitmeplease = findbadfit(time,rxNum,jj);
      if fitmeplease
	
    %if choplet
      %jj
      %time = times{ii}
      %fitmeplease = findbadfit(time,rxNum,jj);
      %if fitmeplease
      %	ww{ii}.AA(:,jj) = findbeamcenter(ww{ii}.x_bin,ww{ii}.y_bin,ww{ii}.map(:,:,jj),cent);   
      %else
	%ww{ii}.AA(:,jj) = NaN(1,7); 
      %end

      % If error with normfit2d, give the struct NaNs
        try
          ww{ii}.AA(:,jj) = ...
	      findbeamcenter(ww{ii}.x_bin,ww{ii}.y_bin,ww{ii}.map(:,:,jj),cent);
        catch err
          ww{ii}.AA(:,jj) = NaN(1,7);
        end % try
      
      else % fitmeplease
        ww{ii}.AA(:,jj) = NaN(1,7);
      end
      
    else % normal demod
      ww{ii}.AA(:,jj) = ...
	  findbeamcenter(ww{ii}.x_bin,ww{ii}.y_bin,ww{ii}.map(:,:,jj),cent);
    end
  end
  
  clear maps x_bin y_bin Afit;
  map = ww{ii}.map;
  x_bin = ww{ii}.x_bin;
  y_bin = ww{ii}.y_bin;
  Afit = ww{ii}.AA;
  
  switch experiment
    case 'bicep2'
      exname = 'b2';
      save(['beammaps/maps_rotated/' char(times{ii}) '_' exname ...
	    '_rot' num2str(dk(ii))],'map','x_bin','y_bin','Afit')
    case 'keck'
      exname = 'keck';  
      if choplet
	save(['beammaps/maps_rotated/' char(times{ii}) '_' exname '_cl_rx' ...
	      num2str(rxNum) '_rot' ...
	      num2str(dk(ii))],'map','x_bin','y_bin','Afit')
      else
	save(['beammaps/maps_rotated/' char(times{ii}) '_' exname '_rx' ...
	      num2str(rxNum) '_rot' ...
	      num2str(dk(ii))],'map','x_bin','y_bin','Afit')
      end
  end

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = findbeamcenter(x_bin,y_bin,map,cent)

winaz = (x_bin > cent.x - 2 & x_bin <= cent.x + 2);
winel = (y_bin > cent.y - 2 & y_bin <= cent.y + 2);
  
x_win = x_bin(winaz);
y_win = y_bin(winel);  

[maskaz maskel] = meshgrid(winaz,winel);
mask = maskaz & maskel;
maptofit = map(winel,winaz);
% No fit if NaNs in map
maptofit(isnan(maptofit)) = 0;

A = normfit2d(x_win,y_win,maptofit,'window');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fitmeplease = findbadfit(time,rxNum,jj)

fitmeplease = 1;
% SUPER STUPID way to get around normfit2d failing on individual maps

if strcmp(time,'20130213_014347')
  if rxNum == 4
    if jj == 275 | jj == 281 | jj == 347 | jj == 363 | jj == 504
      fitmeplease = 0;
    end
  end
end
if strcmp(time,'20130213_092408')
  if rxNum == 4
    if jj == 44 | jj == 45 | jj == 60 | jj == 61 | jj == 77 | jj == 78 ...
	  | jj == 94 | jj == 110 | jj == 111 | jj == 125 | jj == 126 | ...
	  jj == 127 | jj == 223 | jj == 234 | jj == 272 | jj == 280 | ...
	  jj == 327 | jj == 331 | jj == 383 | jj == 467
      fitmeplease = 0;
    end
  end
end
if strcmp(time,'20130219_122110')
  if rxNum == 4
    if jj == 115 | jj == 223 | jj == 275 | jj == 289 | jj == 290 | jj ...
	  == 330 | jj == 364 | jj == 388
      fitmeplease = 0;
    end
  end
end
if strcmp(time,'20130220_211148')
  if rxNum == 4
    if jj == 223
      fitmeplease = 0;
    end
  end
end
if strcmp(time,'20130221_123717')
  if rxNum == 4
    if jj == 1 | jj == 243 | jj == 246 | jj == 278 | jj == 372
      fitmeplease = 0;
    end
  end
end
if strcmp(time,'20130224_031531')
  if rxNum == 4
    if jj == 272 | jj == 275 | jj == 278 | jj == 280 | jj == 380 | jj == ...
	  383 | jj == 391 | jj == 395 | jj == 470 | jj == 504 | jj == 514 
      fitmeplease = 0;
    end
  end
end
if strcmp(time,'20130224_110215')
  if rxNum == 4
    if jj == 238 | jj == 282 | jj == 331 | jj == 372
      fitmeplease = 0;
    end
  end
end
if strcmp(time,'20130301_142525')
  if rxNum == 4
    if jj == 292 | jj == 348 | jj == 366 | jj == 367
      fitmeplease = 0;
    end
  end
end

% 2014
if strcmp(time,'20140301_073616')
  if rxNum == 0
    if jj == 22 | jj == 28 | jj == 38 | jj == 40 | jj == 42 | jj == 44 | jj ...
	  == 45 | jj == 56 | jj == 59 | jj == 62 | jj == 63 | jj == 73 | jj ...
	  == 97 | jj == 107 | jj == 110 | jj == 113 | jj == 159 | jj == 160 ...
	  | jj == 162 | jj == 177 | jj == 179 | jj == 187 | jj == 189 | jj ...
	  == 191 | jj == 193 | jj == 203 | jj == 205 | jj == 206 | jj == ...
	  209 | jj == 221 | jj == 237 | jj == 239 | jj == 245
      fitmeplease = 0;
    end
  end
end

if strcmp(time,'20140301_151619')
  if rxNum == 2
    if jj == 21 | jj == 27 | jj == 43 | jj == 46 | jj == 48 | jj == 54 | jj ...
	  == 80 | jj == 88 | jj == 95 | jj == 97 | jj == 111 | jj == 113 | ...
	  jj == 178 | jj == 195 | jj == 219 | jj == 220 | jj == 245 | jj == ...
	  285 | jj == 286 | jj == 288 | jj == 291 | jj == 294 | jj == 296 | ...
	  jj == 304 | jj == 318 | jj == 319 | jj == 323 | jj == 324 | jj == ...
	  328 | jj == 335 | jj == 346 | jj == 357 | jj == 360 | jj == 361
      fitmeplease = 0;
    end
  end
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%various checks
%{
collmap=A(1).collmap;
dklist=[];
for ii=1:length(collmap)
  dklist=[dklist collmap(ii).dk];
end
%}

%{
imagesc(mapstoadd(ii,jj).map)
hold on
plot(az(ii)/2,el(ii)/2,'xk')
%}
%{

imagesc(x_bin,y_bin,map_same(:,:,ii))
hold on
plot(AA(2,1),AA(3,1),'x')
plot(AA(2,2),AA(3,2),'xy')
plot(AA(2,3),AA(3,3),'xg')
plot(AA(2,4),AA(3,4),'xm')
%}

%{
function [x_bin_padded,y_bin_padded,map_padded]=centermaponbeam(xcent,ycent,x_bin,y_bin,map)
x_res=x_bin(2)-x_bin(1);
numtopad=round((mean(x_bin)-xcent)/x_res)*2-1;
if numtopad>0
  map_padded=padarray(map,[0 numtopad],0,'pre');
  x_bin_padded=x_bin;
  for jj=1:abs(numtopad)
    x_bin_padded=[x_bin_padded(1)-x_res x_bin_padded];
  end
else
  map_padded=padarray(map,[0 abs(numtopad)],0,'post');
  x_bin_padded=x_bin;
  for jj=1:abs(numtopad)
    x_bin_padded=[x_bin_padded x_bin_padded(end)+x_res];
  end
end
y_res=y_bin(2)-y_bin(1);
numtopad=round((mean(y_bin)-ycent)/y_res)*2-1;
if numtopad>0
  map_padded=padarray(map_padded,abs(numtopad),0,'pre');
  y_bin_padded=y_bin;
  for jj=1:abs(numtopad)
    y_bin_padded=[y_bin_padded(1)-y_res y_bin_padded];
  end
else
  map_padded=padarray(map_padded,abs(numtopad),0,'post');
  y_bin_padded=y_bin;
  for jj=1:abs(numtopad)
    y_bin_padded=[y_bin_padded y_bin_padded(end)+y_res];
  end
end
%}

return

%}