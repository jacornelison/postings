function [mapstruct,apsopt,err] = aps_prepare_maps(mapnames,apsopt)
% [mapstruct,apsopt,err] = aps_prepare_maps(mapnames,apsopt)
%
% this is used in reduc_makeaps and reduc_plotcomap_pager to avoid
% duplication
%
% mapnames is a cell - no wild cards!
% the maps are read in and combined as indicated in the apsopt
% apsopt get modified
%
% The expand_ac and scalefac arguments now have two forms - old and
% new:
% In old form expand_ac is non-cell and scalefac is simple
% cell. This is now obsolete but continues to work for the present
% In new form expand_ac is cell over files and scalefac is
% cell-of-cells which first level is over files and second level is
% over components within each cell (noise, lcdm, dust etc)
% The new form is required so that one can mix ac's and map input
% files - see simruncode/make_BK14xWPlanckLT

if ~iscell(mapnames)
  mapnames={mapnames};
end

err=0;
mapstruct=[];

for jj=1:size(mapnames,2)
  try
    disp(['Loading map ',num2str(jj),': ',[apsopt.hostdir{jj},mapnames{jj}]])
    cmapstruct = load([apsopt.hostdir{jj},mapnames{jj}]);
    if ~isfield(cmapstruct,'m')
      disp(['missing field in map structure: ',[apsopt.hostdir{jj},mapnames{jj}]])
      err=1;
      return;
    end
  catch err
    disp(['Failed to load file: ',[apsopt.hostdir{jj},mapnames{jj}]])
    return;
  end

  % expand the m to same size as map
  if(isfield(cmapstruct,'ac'))
    if(~iscell(cmapstruct.ac))
      cmapstruct.m=repmat(cmapstruct.m,[size(cmapstruct.ac,1),1]);
    else
      cmapstruct.m=repmat(cmapstruct.m,[size(cmapstruct.ac{1},1),1]);
    end
  else
    cmapstruct.m=repmat(cmapstruct.m,[size(cmapstruct.map,1),1]);
  end

  % strip to first half of jack split
  if jj>1
    if(isfield(apsopt,'firsthalf'))
      cmapstruct.ac=cmapstruct.ac(:,1)
    end
    if(isfield(apsopt,'secondhalf'))
      cmapstruct.ac=cmapstruct.ac(:,2)
    end
  end
  
  % Calibrate here before coadding over rx/making map
  if length(apsopt.ukpervolt)==1
    cmapstruct=cal_maps(cmapstruct,apsopt.ukpervolt{1},apsopt);
  else
    cmapstruct=cal_maps(cmapstruct,apsopt.ukpervolt{jj},apsopt);
  end

  % if requested coadd files over receivers - else leave them in rx format
  if isfield(cmapstruct, 'ac')
    % this is confusing, because apsopt.coaddrx is badly named, change to
    % apsopt.docoadd instead? Currently coaddrx=0 still makes a coadd to
    % per-rx size - should do nothing instead, i.e. docoadd=0
    switch apsopt.coaddrx
      case 2
        cmapstruct.ac = coadd_ac_overfreq(cmapstruct.ac,cmapstruct.coaddopt);
      case 1
        cmapstruct.ac = coadd_ac_overrx(cmapstruct.ac);
      case 0
        cmapstruct.ac = coadd_ac_intorx(cmapstruct.ac,cmapstruct.coaddopt);
      otherwise
        disp('Do no additional coaddition')
    end
    % Assume that all m structures are equivalent since ac have been
    % coadded together if length(ac) ~= length(m)
    if iscell(cmapstruct.ac)
      nac = size(cmapstruct.ac{1}, 1);
    else
      nac = size(cmapstruct.ac, 1);
    end
    if nac ~= length(cmapstruct.m)
      cmapstruct.m = cmapstruct.m(1:nac,1);
    end
    clear nac
  end
  
  if strcmp(apsopt.howtojack,'none')
    % don't jackknife - we want seperate spectra for each jack split
    % the jackknife dimension is also used for jack0 stuff (abscal, deltaaps)
    % so if we want cross spectra between the two halfs of a split, these
    % all need to go along the first dimension of the map structure. So
    % flatten map out: 
    % J1_100 J2_100
    % J1_150 J2_150
    % --> J1_100; J2_100; J1_150; J2_150
    % etc.
    % if coaddopt.jacktype~='0'
    if isfield(cmapstruct, 'ac')
      if iscell(cmapstruct.ac)
        for i = 1:numel(cmapstruct.ac)
          cmapstruct.ac{i} = flatten_ac(cmapstruct.ac{i});
        end
      else
        cmapstruct.ac = flatten_ac(cmapstruct.ac);
      end
    else
      if iscell(cmapstruct.map)
        for i = 1:numel(cmapstruct.map)
          cmapstruct.map{i} = flatten_ac(cmapstruct.map{i});
        end
      else
        cmapstruct.map = flatten_ac(cmapstruct.map);
      end
    end
  end
  
  % select/use multiple times the observations (freqs) indicated 
  % with the expand_ac field
  % NEW VERSION: expand_ac is now a cell over files so can have
  % some files be ac's and some maps and still have this work -
  % once this is common practice delete the obsolete version below
  if(isfield(apsopt,'expand_ac'))
    if(iscell(apsopt.expand_ac))
      if(~isempty(apsopt.expand_ac{jj}))
        disp(sprintf('expanding ac from this file as requested'));
        if(iscell(cmapstruct.ac))
          for ii = 1:numel(cmapstruct.ac)
            cmapstruct.ac{ii} = cmapstruct.ac{ii}(apsopt.expand_ac{jj});
          end
        else
          cmapstruct.ac = cmapstruct.ac(apsopt.expand_ac{jj});
        end
        % expand the m as well
        cmapstruct.m=cmapstruct.m(apsopt.expand_ac{jj});
      end
    end
  end
  
  % if requested scale component ac's - for dust or noise scaling
  % etc
  % NEW VERSION: scalefac is now a cell over files of cell arrays
  % so can have some files be ac's and some maps and still have
  % this work - once this is common practice delete the obsolete
  % version below.
  % Note that currently scaling of maps is not implemented
  if(isfield(apsopt,'scalefac'))
    if(iscell(apsopt.scalefac{1}))

      % if we are doing scaling all ac's become cell
      if(isfield(cmapstruct,'ac'))
        if(~iscell(cmapstruct.ac))
          cmapstruct.ac = {cmapstruct.ac};
        end
      end

      if(~isempty(apsopt.scalefac{jj}))
      
        disp('applying scaling to ac from this file as requested');
        % loop over each cell component of ac (noise,signal etc)
        % set cal_coadd_ac to only scale the signal portion, not the variance
    
        for i=1:numel(cmapstruct.ac)
          % zero results in blank maps - this doesn't
          apsopt.scalefac{jj}{i}(apsopt.scalefac{jj}{i}==0)=1e-99;
          % apply separate scale factor for each band
          for j=1:size(cmapstruct.ac{i},1)
            % and to both halfs of a jack if applies
            for k=1:size(cmapstruct.ac{i},2)
              cmapstruct.ac{i}(j,k)=cal_coadd_ac(cmapstruct.ac{i}(j,k),apsopt.scalefac{jj}{i}(j),cmapstruct.coaddopt,0,1);
            end
          end
        end
        
      end
    end
  end
  
  % no merging of the structures when on the first file, 
  % unless it is the only file (auto spectrum)
  if jj==1
    mapstruct = cmapstruct;
    if isfield(mapstruct,'ac')
      if(~iscell(mapstruct.ac))
        nmap(jj)=size(mapstruct.ac,1);
      else
        nmap(jj)=size(mapstruct.ac{1},1);
      end
    else
      nmap(jj)=size(mapstruct.map,1);
    end
    if length(mapnames)>1
      continue;
    end
  end
  
  % the next if statement just applies if we read more than one file
  if length(mapnames)>1
    
    % keep all the m's
    mapstruct.m=struct_merge(mapstruct.m,cmapstruct.m);
    
    % normally files contain ac's - if one or more contains
    % maps we convert all to maps within loop over files - this
    % should probably be made the default but there is complex
    % stuff acting on ac's below which will break if we do this...
    %
    % require mapstruct to not already have map; this circumvent cases
    % where mapstruct.map exists which come from having merged with 
    % cmapstruct that does not have ac structure but only has map.
    % So mapstruct.ac has N items, where mapstruct.map has N+1 items.
    % If the next (N+2th) cmapstruct only has ac and is merged with mapstruct.ac,
    % the merged mapstruct would now have N+1 ac's, but N+2 m's. And the 
    % N+1th ac would get the N+1th map slot instead of N+2th map slot where
    % it should be.
    % may break cases ~iscell(campstruct.ac); beware.
    if (isfield(cmapstruct,'ac')&isfield(mapstruct,'ac')&~isfield(mapstruct,'map'))
      % this handles the special case where cross spectra between different maps
      % types (for instance S+N & B) are calculated. In this case one ac structure
      % may be a cell array the other one not
      if iscell(cmapstruct.ac)~=iscell(mapstruct.ac)
        if ~isfield(mapstruct,'map')
          mapstruct.map=make_map(mapstruct.ac,mapstruct.m,mapstruct.coaddopt);
        end
        cmapstruct.map =make_map(cmapstruct.ac,cmapstruct.m,cmapstruct.coaddopt);
        nmap(jj)=size(cmapstruct.map,1);
        mapstruct.map = struct_merge(mapstruct.map,cmapstruct.map);
      else 
        % concatenate the ac structures over files
        if(~iscell(cmapstruct.ac))
          nmap(jj)=size(cmapstruct.ac,1);
          % the struct_merge will fail for cross half jacks where
          % one map is full the other one is jack. In this case
          % fill the missing bit of the structure with an empty
          % ac. This will cause the jackknife map part to skip
          % over it and leave the non-jack map in place
          % will do nothing when the second dimensions are the same.
          % this enables Full map x Jack map
          [mapstruct.ac,cmapstruct.ac]=equalize_ac_2nd_dim(mapstruct.ac,cmapstruct.ac);
          
          mapstruct.ac=struct_merge(mapstruct.ac,cmapstruct.ac);
        else
          nmap(jj)=size(cmapstruct.ac{1},1);
          for i=1:length(cmapstruct.ac)
            % see comment above
            [mapstruct.ac{i},cmapstruct.ac{i}]=equalize_ac_2nd_dim(mapstruct.ac{i},cmapstruct.ac{i});
            
            mapstruct.ac{i}=struct_merge(mapstruct.ac{i},cmapstruct.ac{i});
          end
        end
      end
    else % One or more files contains maps rather than ac's
      % make sure we got maps
      if ~isfield(mapstruct,'map')
        mapstruct.map=make_map(mapstruct.ac,mapstruct.m,mapstruct.coaddopt);
      end
      if ~isfield(cmapstruct,'map') 
        cmapstruct.map=make_map(cmapstruct.ac,cmapstruct.m,cmapstruct.coaddopt);
      end
      nmap(jj)=size(cmapstruct.map,1);
      % concatenate the map structures
      % if B-mode map exists in only one of cmapstruct or mapstruct, keep it.
      % default behavior of struct_merge keeps common fields only.
      if isfield(mapstruct.map,'B')&~isfield(cmapstruct.map,'B')
         cmapstruct.map.B=[];
      end
      if isfield(cmapstruct.map,'B')&~isfield(mapstruct.map,'B')
        tmp = mapstruct.map;
        for xx=1:size(mapstruct.map,1)
          tmp(xx).B=[];
        end
        mapstruct.map = tmp;
        clearvars tmp
      end 
      mapstruct.map = struct_merge(mapstruct.map,cmapstruct.map);
    end
  end
  
end % loop over the maps to be crossed

% if requested select a component (signal or noise etc.)
% from a combined map. (this is a dirty trick; in the mainline
% sims this should never be necessary)
if(isfield(apsopt,'select_ac'))
  if(iscell(mapstruct.ac))
    mapstruct.ac=mapstruct.ac{apsopt.select_ac};
  end
end

% select/use multiple times the observations (freqs) indicated 
% with the expand_ac field
% OBSOLETE: make expand_ac a cell over files so can have some files
% be ac's and some maps and still have this work - see inside loop
% over files above
if(isfield(apsopt,'expand_ac'))
  if(~iscell(apsopt.expand_ac))
    if(iscell(mapstruct.ac))
      for ii = 1:numel(mapstruct.ac)
        mapstruct.ac{ii} = mapstruct.ac{ii}(apsopt.expand_ac);
      end
    else
      mapstruct.ac = mapstruct.ac(apsopt.expand_ac);
    end
    % expand the m as well
    mapstruct.m=mapstruct.m(apsopt.expand_ac);
  end
end

% if requested scale component maps before they are added
% together in make_map
% OBSOLETE: scalfac should now be a cell-of-cells so can have some
% files be ac's and some maps and still have this work - see inside
% loop over files above
if(isfield(apsopt,'scalefac'))
  if(~iscell(apsopt.scalefac{1}))
    if(~iscell(mapstruct.ac))
      mapstruct.ac = {mapstruct.ac};
    end
    
    disp('applying scaling to maps');
    % loop over each cell component of ac (noise,signal etc)
    % set cal_coadd_ac to only scale the signal portion, not the variance
    
    for i=1:numel(mapstruct.ac)
      % zero results in blank maps - this doesn't
      apsopt.scalefac{i}(apsopt.scalefac{i}==0)=1e-99;
      % apply separate scale factor for each receiver
      for j=1:size(mapstruct.ac{i},1)
        % and to both halfs of a jack if applies:
        for k=1:size(mapstruct.ac{i},2)
          mapstruct.ac{i}(j,k)=cal_coadd_ac(mapstruct.ac{i}(j,k),apsopt.scalefac{i}(j),mapstruct.coaddopt,0,1);
        end
      end
    end
  end
end

% coadd over multiple receivers if requested
if (iscell(apsopt.overall) || apsopt.overall) && isfield(mapstruct, 'ac') && ~isempty(apsopt.overall)
  if iscell(apsopt.overall)
    if iscell(mapstruct.ac)
      %loop over each cell component of ac, that is the types
      for jj=1:numel(mapstruct.ac)
        clear actmp
        for kk=1:length(apsopt.overall)
          actmp(kk,:)=coadd_ac_overrx(mapstruct.ac{jj}(apsopt.overall{kk},:));
        end
        mapstruct.ac{jj}=actmp;
      end
    else %single map
      clear actmp
      for kk=1:length(apsopt.overall)
        actmp(kk,:)=coadd_ac_overrx(mapstruct.ac(apsopt.overall{kk},:));
      end
      mapstruct.ac=actmp;
    end
  else % combining everything in this case.  coadd_ac_overrx can handle cells
    mapstruct.ac=coadd_ac_overrx(mapstruct.ac);
  end
end    

% boost sky coverage so keck maps can match B3
if(isfield(apsopt,'match_sky_coverage'))
  disp('changing sky coverage (map defn) as instructed');
  m1idx = rvec(apsopt.match_sky_coverage{1});
  m1a = mapstruct.m(m1idx);
  if isnumeric(apsopt.match_sky_coverage{2})
    m2a = mapstruct.m(apsopt.match_sky_coverage{2});
  else
    m2a = apsopt.match_sky_coverage{2};
  end
  if numel(m2a) == 1 && numel(m2a) < numel(m1a)
    m2a = repmat(m2a, 1, numel(m1a));
  end
  if iscell(mapstruct.ac)
    nac2 = size(mapstruct.ac{1}, 2);
  else
    nac2 = size(mapstruct.ac, 2);
  end
  for i=1:numel(m1idx)
    m1=m1a(i);
    m2=m2a(i);
    for j=1:nac2
      if iscell(mapstruct.ac)
        for k=1:length(mapstruct.ac)
          mapstruct.ac{k}(m1idx(i),j) = expandac(m1, m2, ...
              mapstruct.ac{k}(m1idx(i),j));
        end
      else
        mapstruct.ac(m1idx(i),j) = expandac(m1, m2, mapstruct.ac(m1idx(i),j));
      end
    end
    mapstruct.m(m1idx(i))=m2;
  end
end

% convert from ac to map already
if ~isfield(mapstruct,'map')
  mapstruct.map=make_map(mapstruct.ac,mapstruct.m,mapstruct.coaddopt);
end

% expand polrot appropriately
if(isfield(apsopt,'polrot') & ~isempty(apsopt.polrot))
  switch length(apsopt.polrot)
   case 1
    % single value specified
    apsopt.polrot=repmat(apsopt.polrot,1,sum(nmap));
   case size(mapstruct.map,1)
    % one value per map - don't need to do anything
    if (iscell(apsopt.overall) || apsopt.overall)
      warning(['Polrot assumed per map after application of apsopt.overall.'])
    end
   case size(mapnames,2)
    % one value per file
    if (iscell(apsopt.overall) || apsopt.overall)
      warning(['The application of polrot does not automatically handel the tweak to the maps done with apsopt.overall.'])
    end
    polrot=apsopt.polrot;
    apsopt.polrot=repmat(apsopt.polrot(1),1,nmap(1));
    for jj=2:length(mapnames)
      apsopt.polrot=[apsopt.polrot,repmat(polrot(jj),1,nmap(jj))];
    end
   otherwise
    error('dont know how to expand apsopt.polrot')
  end
end

% apply polrot
if(isfield(apsopt,'polrot') & ~isempty(apsopt.polrot))
  mapstruct.map=rotqumaps(mapstruct.map,apsopt.polrot);
end

% this kludge is necessary to get rid of cell coaddopt
if iscell(mapstruct.coaddopt)
  mapstruct.coaddopt=mapstruct.coaddopt{1};
else
  mapstruct.coaddopt=mapstruct.coaddopt;
end

% de-rotate Q and U if the projection is anything other than RA/DEC
switch mapstruct.coaddopt.coaddtype
case {0,1,2,5}
  for i=1:size(mapstruct.map,1)
    for j=1:size(mapstruct.map,2)
      if ~isempty(mapstruct.map(i,j).Q)
        [mapstruct.map(i,j).Q,mapstruct.map(i,j).U]=derotate_qu(mapstruct.map(i,j).Q,mapstruct.map(i,j).U,mapstruct.m(i));
      end
    end
  end
end

% jackknife the maps
switch apsopt.howtojack
 case 'dim2'
  % normal case - jackknife along 2nd map dim
  mapstruct.map=jackknife_map(mapstruct.map);

 case 'dim1'
  % take diff along the 1st map dim
  % for rx exp jackknives, we need to call jackknife_map multiple times
  % the first part is here a dummy place holder for where usually
  % the auto spectra are, these maps will be zero
  for i=1:size(mapstruct.map,1)
    mapi=[mapstruct.map(i,1),mapstruct.map(i,1)]';
    dmap(i,1)=jackknife_map(mapi,1);
    dm(i,1)=mapstruct.m(i);
  end
  c = size(mapstruct.map,1)+1;
  % in the same way cross spectra are done also the subtraction
  % is set up:
  for i=1:size(mapstruct.map,1)-1
    for j=i+1:size(mapstruct.map,1)
      mapi=[mapstruct.map(i,1),mapstruct.map(j,1)]';
      dmap(c,1)=jackknife_map(mapi,1);  
      dm(c,1)=mapstruct.m(j);
      c=c+1;
    end
  end
  mapstruct.map=dmap;
  mapstruct.m=dm;
  % keep a record that this is jack spectrum
  mapstruct.coaddopt.jacktype='f';
end

% Calculate the apodization masks
% note: Qvar and Uvar are not now passed forward smoothed which is
% a change from previous behavior
mapstruct.map=add_masks(mapstruct.m(1),mapstruct.map,apsopt.smoothvarmaps);

% special option to force in a mask from a different file
if(isfield(apsopt,'alt_mask_filename'))
  % load alternate mask from given filename
  load(apsopt.alt_mask_filename);
  disp(sprintf('forcing in mask from file %s',apsopt.alt_mask_filename));
  if(~all(size(map.Tw)==size(mapstruct.map(1).Tw)))
    if(size(map.Tw,1)>size(mapstruct.map(1).Tw,1))
      % attempt to downsample the alt mask
      disp('downsampling the alt mask to match map(1)');
      dsf=size(map.Tw,1)./size(mapstruct.map(1).Tw,1);
      map.Tw=map.Tw(1:dsf:end,1:dsf:end);
      map.Pw=map.Pw(1:dsf:end,1:dsf:end);
    else
      % attempt to upsample the alt mask
      disp('upsampling the alt mask to match map(1)');
      [x,y]=meshgrid(m.x_tic,m.y_tic);
      map.Tw=interp2(x,y,map.Tw,mapstruct.m.x_tic,mapstruct.m.y_tic');
      map.Pw=interp2(x,y,map.Pw,mapstruct.m.x_tic,mapstruct.m.y_tic');
    end
  end
  % copy the alt mask to all maps
  disp(sprintf('copying this alt mask to all %d fields of map',numel(mapstruct.map)));
  for i=1:numel(mapstruct.map)
    mapstruct.map(i).Tw=map.Tw;
    mapstruct.map(i).Pw=map.Pw;    
  end
end

% saw off the mask of some maps to allow experiments with
% non-overlapping sky coverage
if(isfield(apsopt,'saw_off_mask'))
  disp('sawing off masks for some bands as instructed');
  x=mapstruct.m(1).x_tic;
  % these numbers say where roll-up starts and ends
  [dum,i1]=min(abs(x-apsopt.saw_off_mask{1}(1)));
  [dum,i2]=min(abs(x-apsopt.saw_off_mask{1}(2)));
  s(1:i1)=0;
  s(i1:i2)=sind(linspace(0,90,i2-i1+1)).^2;
  s(i2:length(x))=1;
  m=ones(mapstruct.m(1).ny,1)*s;
  for i=apsopt.saw_off_mask{2}
    for j=1:size(mapstruct.map,2)
      mapstruct.map(i,j).Tw=mapstruct.map(i,j).Tw.*m;
      mapstruct.map(i,j).Pw=mapstruct.map(i,j).Pw.*m;
    end
  end
end

% if requested construct and/or copy common mask
if apsopt.commonmask
  mapstruct.map=insert_common_mask(mapstruct.m,mapstruct.map,apsopt.commonmask,apsopt.commonmasksel);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mapstruct=cal_maps(mapstruct,ukpervolt,apsopt)

% cal the maps unless they are already in uK.
if strcmp(ukpervolt, 'none')
  warning('do not apply the abs cal here')
  return
end

% presumably always will have a codaddopt...
if(isfield(mapstruct,'coaddopt'))

  if iscell(mapstruct.coaddopt)
    coaddopt=mapstruct.coaddopt{1};
  else
    coaddopt=mapstruct.coaddopt;
  end
  
  if(isfield(coaddopt,'ukpv_applied')&&~isempty(coaddopt.ukpv_applied))
    disp('The map was already calibrated during coadd.')
    ukpervolt=ones(size(ukpervolt));
    
    % find out if we are a sim
  else
    if(isfield(coaddopt.mapopt{1},'simopt'))

      % fetch the ukpervolt used to make the sim
      if ~iscell(coaddopt.mapopt{1}.simopt)
        sukpervolt=coaddopt.mapopt{1}.simopt.ukpervolt;
      else
        % not sure that this ever happens...
        sukpervolt=coaddopt.mapopt{1}.simopt{1}.ukpervolt;
      end

      % if we are a sim and ac is composite (e.g. signal&noise) then
      % the situation is dangerous...
      if(isfield(mapstruct,'ac'))
        if(iscell(mapstruct.ac))
          if(ukpervolt~=sukpervolt)
            % if the ukpervolt specified in simopt does not match the
            % apsopt value then we don't really know what to do - we
            % don't really know which of the cell elements of ac should use
            % the simopt.ukpervolt value and which the
            % apsopt.ukpervolt value.
            error('ac is cell and the ukpervolt in simopt is not the same as the one in apsopt - the information is not available to handle this');
          end
        end
      end 
      % fetch the ukpervolt used to make this sim
      if ~apsopt.overrideukpervolt
        warning('ukpervolt forced to that used to make this sim');
        ukpervolt=sukpervolt
      else
	warning('overriding sim ukpervolt replacement, retaining apsopt ukpervolt')
      end
    end
  end
else
  warning('No coaddopt found in map structure.')
end

% apply the ukpervolt
if ~isfield(mapstruct,'map')
  if(~isfield(apsopt,'multac'))
    mapstruct.ac=cal_coadd_ac(mapstruct.ac,ukpervolt,mapstruct.coaddopt);
  else
    % This is special switch for making normalized coverage sims
    % We have made a coadd of "subsim" style weight=1, var=1
    % pairmaps. We wish to scale out any ukpervolt used for the
    % signal portion (could have been 1 but may not have been) and
    % then multiply the entire ac structure by some empirical scale
    % factor. This scale factor is adjusted to match the var maps
    % of a standard real coadd, and can then be used in projection
    % sims.
    disp('Special apsopt.multac switch in effect - will cal only signal portion of ac and then mult whole ac by requested factors');
    mapstruct.ac=cal_coadd_ac(mapstruct.ac,ukpervolt,mapstruct.coaddopt,0,1);
    mapstruct.ac=multac(mapstruct.ac,apsopt.multac);
  end
else
  mapstruct.map=cal_coadd_maps(mapstruct.map,ukpervolt);
end

return

function [ac1,ac2]=equalize_ac_2nd_dim(ac1,ac2)
% function [ac1,ac2]=equalize_ac_2nd_dim(ac1,ac2)
% have both ac structures be size(ac,2)=2 if either
% ac1 or ac2 has size(ac,2)=2. The missing bits are
% filled with empty fields
if size(ac1,2)==2 && size(ac2,2)==1
  for i=1:size(ac2,1)
    ac2(i,2)=get_empty_ac(ac2(i,1));
  end
elseif size(ac1,2)==1 && size(ac2,2)==2
  for i=1:size(ac1,1)
    ac1(i,2)=get_empty_ac(ac1(i,1));
  end
end
return

function ac_flat = flatten_ac(ac)
% function ac_flat = flatten_ac(ac)
% ac not a cell here!
% size(ac) = 2x2 --> 4x1
if size(ac,2)==2
  ac_t = ac';
  for i=1:numel(ac_t)
    cac(i)=ac_t(i);
  end
  ac_flat = cac';    
else
  ac_flat = ac;
end
return
