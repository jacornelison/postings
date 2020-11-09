function ffbm_resizemaps(experiment,year,mapname,concat,cutmaps,addnoise,renorm,makeinnermap,relgainmap,savename)
% ffbm_resizemaps(experiment,year,mapname,concat,cutmaps,addnoise,renorm,makeinnermap,relgainmap,savename)
%
% Resize, concatenate, add noise, renormalize maps
% Make 1.2 deg radius beam maps, 
% add relgain into beam maps
% Unlikely to just work
%
% experiment = 'bicep2','keck'
%
% mapname = string of map to load
%           should contain matrix map, struct ad
%           can be left empty
%
% concat = 1 to concatenate keck composite maps
%          takes the output of compositemaps.m
%          concatenates the maps into one file
%
% cutmaps = size to cut map into
%         = 0 to not cut the maps
%
% addnoise = 0 to not add noise into the map
%          = 1 to add white noise into beam map
%             the noise added is at the level of each individual map
%             as estimated from estimatenoise.m
%          = 2 to add annulus mean and noise into map
%
% makeinnermap = 0 to not do anything
%              = 1 to make inner 1.2 deg map. uses a mean filter
%                and zeros things outside r=1.2 degs
%
% renorm = 1 to integral renormalize map
%        = 0 to not renormalize map
%
% relgainmap = 0 to not do anything
%            = 2012 or 2013 to insert relgain. Should probably look at it
%              carefully, it's unlikely to just work
%
% savename = filename to save
%
% resizemaps('keck','maps/keck_8deg',1,0,1,0,0,'maps_wnoise/keck_8deg')
% resizemaps('keck','maps/rx1_all',0,0,1,0,0,'maps_wnoise/rx1_all')
% resizemaps('keck','maps/rx1_all',0,6.1,1,0,0,'maps_wnoise/rx1_all_6deg')
% resizemaps('keck',[],1,0,0,0,0,'maps/keck_8deg')
% resizemaps('keck',[],1,0,1,1,0,0,'maps_wnoise/keck_8deg')
% resizemaps('keck',[],1,6.1,1,1,0,'maps_wnoise/keck_6deg')
% resizemaps('keck','maps_wnoise/keck_6deg',0,3.1,0,1,0,'maps_wnoise/keck_3deg')
% resizemaps('keck','maps/keck_8deg',0,0,2,0,0,'keckmaps_wnoise/keck2012_8deg_wnoise_pknorm_2')
% resizemaps('keck','keckmaps_wnoise/keck2012_8deg_wnoise_pknorm',0,0,0,1,0,'keckmaps_wnoise/keck2012_8deg_wnoise')
% resizemaps('keck','keckmaps_wnoise/keck2012_8deg_wnoise',0,0,0,0,1,'keckmaps_wnoise/keck2012_8deginner_wnoise')
% resizemaps('keck','maps/keck_8deg',0,0,2,1,0,'keckmaps_wnoise/keck2012_8deg_wnoise_med')
% resizemaps('keck','maps/keck_8deg',0,0,2,0,0,'keckmaps_wnoise/keck2012_wnoise_pknorm')
% resizemaps('keck','keckmaps_wnoise/keck2012_wnoise_pknorm',0,0,0,1,0,'keckmaps_wnoise/keck2012_wnoise')
% ffbm_resizemaps('keck',[],1,0,0,0,0,'maps/keck2013_all')
% ffbm_resizemaps('keck','maps/keck2013_all',0,0,2,0,0,'keckmaps_wnoise/keck2013_wnoise_pknorm')
% ffbm_resizemaps('keck','keckmaps_wnoise/keck2013_wnoise_pknorm',0,0,0,1,0,'keckmaps_wnoise/keck2013_wnoise')
% ffbm_resizemaps('keck','keckmaps_wnoise/keck2013_wnoise',0,0,0,0,1,'keckmaps_wnoise/keck2013_wnoise_inner')
% ffbm_resizemaps('keck','keckmaps_wnoise/keck2013_wnoise',0,0,0,0,0,2013,'keckmaps_wnoise/keck2013_wnoise_rg')
% ffbm_resizemaps('keck','keckmaps_wnoise/keck2012_wnoise',0,0,0,0,0,2012,'keckmaps_wnoise/keck2012_wnoise_rg')
% ffbm_resizemaps('keck','keckmaps_wnoise/keck2012_wnoise_blnorm',0,0,0,0,1,0,'keckmaps_wnoise/keck2012_wnoise_blnorm_inner')

if ~exist('addnoise','var')
  addnoise = 0;
  renorm = 0;
end
if isempty(mapname)
  mapname = 0;
end

if ischar(mapname)
  load(mapname)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concatenate maps 
if concat
  for ii = 1:5
    rxNum = ii - 1;
    a{ii} = load(['beammaps/maps_cropped/rx' num2str(rxNum) '_' ...
	  num2str(year) '_all.mat']);
  end
  
  map = cat(3,a{1}.map,a{2}.map,a{3}.map,a{4}.map,a{5}.map);
  ad = a{1}.ad;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Crop maps
if cutmaps
  sizeindeg = cutmaps;
  x_bin = ad.t_val_deg{1};
  y_bin = ad.t_val_deg{2};
  winaz = x_bin >= -sizeindeg/2 & x_bin <= sizeindeg/2;
  winel = y_bin >= -sizeindeg/2 & y_bin <= sizeindeg/2;

  win_x = find(winaz);
  win_x_bin = x_bin(win_x);
  win_y = find(winel);
  win_y_bin = y_bin(win_y);

  map_new = map(win_y(1):win_y(end),win_x(1):win_x(end),:);

  clear map
  map = map_new;
  
  % Make ad structure
  sizeindeg_x = length(win_x)*0.1; %(win_x_bin(2)-win_x_bin(1))
  sizeindeg_y = length(win_y)*0.1; %(win_y_bin(2)-win_y_bin(1))
  Field_size_deg = [sizeindeg_x sizeindeg_y ];
  N_pix = [length(win_x) length(win_y)];
  ad = calc_ad2(Field_size_deg,N_pix);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding noise to map
if addnoise == 1
  mapfilename.map = map;
  mapfilename.x_bin = ad.t_val_deg{1};
  mapfilename.y_bin = -ad.t_val_deg{2};
  noise = ffbm_estimatenoise(mapfilename,'keck',2012,'all',0,0);
  
  map_new = NaN(size(map));
  for ii = 1:size(map,3)
    maptmp = reshape(map(:,:,ii),1,[]);
    % Generate a different noise level for each detector
    noisemap = normrnd(0,noise(ii),size(maptmp));
    nanindex = find(isnan(maptmp));
    maptmp(nanindex) = noisemap(nanindex);
    map_new(:,:,ii) = reshape(maptmp,size(map,1),size(map,2));
  end
  
  clear map;
  map = map_new;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Add noise + mean of values at that annulus
elseif addnoise == 2
  mapfilename.map = map;
  mapfilename.x_bin = ad.t_val_deg{1};
  mapfilename.y_bin = -ad.t_val_deg{2};
  % noise = ffbm_estimatenoise(mapfilename,'keck',2012,'all',0,0);
  % figure out what radius things are at.
  [xx,yy] = meshgrid(mapfilename.x_bin,mapfilename.y_bin);
  rr = sqrt(xx.^2 + yy.^2);
  delta = 0.1; % half deg steps to start
  numstep = round(max(max(rr))/delta);
  r1 = reshape(rr,1,[]);
  posval = NaN(size(r1));
  for ii = 1:numstep
    %ii
    %keyboard
    binlim = [ii-1 ii]*delta;
    winindex = find(r1 >= binlim(1) & r1 < binlim(2));
    posval(winindex) = ii;
  end
  posvalmap = reshape(posval,length(mapfilename.y_bin),[]);
  
  % Calculate the mean and stddev
  s = NaN(numstep,size(map,3));
  m = NaN(numstep,size(map,3));
  mm = NaN(numstep,size(map,3));
  for jj = 1:size(map,3)
    maptmp = reshape(map(:,:,jj),1,[]);
    for ii = 1:numstep
      winindex = find(posval == ii);
      binvalue = maptmp(winindex);
      binvalue(isnan(binvalue)) = [];
      m(ii,jj) = median(binvalue);
      mm(ii,jj) = mean(binvalue);
      %s(ii,jj) = noise(jj);
      s(ii,jj) = std(binvalue);      
    end
  end
  smedian = nanmedian(s(21:end,:),1);
  s = repmat(smedian,size(s,1),[]);
  
  meanmap = NaN(size(reshape(map,1,[],size(map,3))));
  stdmap = NaN(size(meanmap));
  for ii = 1:numstep
    winindex = find(posval == ii);
    tmp = reshape(m(ii,:),1,1,[]);
    meanmap(1,winindex,:) = repmat(tmp,[1,size(winindex,2),1]);
    tmp = reshape(s(ii,:),1,1,[]);
    stdmap(1,winindex,:) = repmat(tmp,[1,size(winindex,2),1]);;
  end
  %  meanmapmap = reshape(meanmap,size(map,1),[],size(map,3));
  %  stdmapmap = reshape(stdmap,size(map,1),[],size(map,3));
  noisemap = normrnd(meanmap,stdmap);
  
  % Add value into blank spaces in map
  % Use the mean and std from each annulus
  map_new = NaN(size(map));
  for jj=1 : size(map,3)
    maptmp = reshape(map(:,:,jj),1,[]);
    nanindex = find(isnan(maptmp));
    maptmp(nanindex) = noisemap(1,nanindex,jj);
    map_new(:,:,jj) = reshape(maptmp,size(map,1),size(map,2));
  end
  
  clear map;
  map = map_new;
end

%%%%%%%%%%%%%%%%%%%%%%%
if makeinnermap
  % Figure out what the radius is for each pixel.
  x_bin = ad.t_val_deg{1};
  y_bin = ad.t_val_deg{2};
  [xx,yy] = meshgrid(x_bin,y_bin);
  rr = sqrt(xx.^2 + yy.^2);
  
  maptmp = reshape(map,1,[],size(map,3));
  r1 = reshape(rr,1,[]);
  %zeroindex = find(r1 > 1.2);
  zeroindex = find(r1 > 2);
  maptmp(1,zeroindex,:) = 0;
  
  map_new = reshape(maptmp,size(map,1),[],size(map,3));
  clear map
  map = map_new;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Renormalize map
if renorm
  
  map_new = NaN(size(map));
  % Normalize map
  for ii = 1:size(map,3)
    map_new(:,:,ii) = map(:,:,ii)/sum(sum(map(:,:,ii)));
  end
  clear map
  map = map_new;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relgain map
if relgainmap
 
  if strcmp(experiment,'keck')
    if relgainmap == 2012
      [p,ind] = get_array_info('20120303',[],[],[],[],[],[],'obs')
    elseif relgainmap == 2013
      [p,ind] = get_array_info('20130303',[],[],[],[],[],[],'obs')
    end
  end
  scaleA = p.ukpv(ind.b)./p.ukpv(ind.a); 
  scaleA(isnan(scaleA)) = 1;
  
  map_new = map;
  scaleA = reshape(scaleA,[1,1,length(ind.a)]);
  scaleA = repmat(scaleA,[length(ad.t_val{2}),length(ad.t_val{1}),1]);
  map_new(:,:,ind.a) = map(:,:,ind.a).*scaleA;
  map_new(:,:,ind.b) = map(:,:,ind.b);
 
  clear map;
  map = map_new;
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put map in correct structure for reduc_makesim
for ii = 1:size(map,3)
  tmpmap(ii).T = map(:,:,ii);
  % Eventually put Q/U beams here too
end
map = tmpmap;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save

save(['beammaps/composites/' savename],'map','ad');

return
