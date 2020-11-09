function map=get_sig_map(simopt,rlz)
% map=get_sig_map(simopt,rlz)
%
% function for cutting map down to needed size
% initially a subfunction of reduc_makesim.m

disp('generating or loading the signal map...')

if ~iscell(simopt.sigmapfilename) & strcmp(simopt.sigmaptype,'healmap')
  simopt.sigmapfilename={simopt.sigmapfilename};
end

if rlz
  for i=1:numel(simopt.sigmapfilename)
    fn=simopt.sigmapfilename{i};
    indx=strfind(fn,'xxxx');
    if(~isempty(indx))
      simopt.sigmapfilename{i}=sprintf('%s%04d%s',fn(1:indx-1),rlz,fn(indx+4:end));
    else
      error('simopt.sigmapfilename not fed in with proper xxxx format for nrlz>0')
    end
  end
end

switch simopt.sigmaptype

  case 'corrmap'
    % flat sky 2D maps
  
    % make the corrmap filename
    simfile=sprintf('data/%s/corrmap%06d',simopt.sernum,n);  
    simopt.simfile=simfile;
  
    % check if already exists
    d=dir(sprintf('data/%s',simopt.sernum));
    if(~any(strcmp({d.name},sprintf('corrmap%06d.mat',n))))
      % if it doesn't generate it
      reduc_makecorrmap(n,simopt.type,simopt.sernum,simopt.corrmap_apsfile);
    end
  
    % load the sim map
    load(simfile);
    simopt.sigmapfilename=simfile;
    if(~strcmp(simopt.corrmap_apsfile,map.apsfile))
      error(['simopt.corrmap_apsfile doesnt match what is actually in the ' ...
        'corrmap file']) 
    end
  
    % find center of field      
    m=get_map_defn(simopt.type);
    map.racen=mean(m.x_tic); map.deccen=mean(m.y_tic);
  
    % convert from flat sky map to ra/dec map
    % - this is an approximation only valid for small maps!
    map.y_tic=map.ad.t_val_deg+map.deccen;
    map.x_tic=map.racen+map.ad.t_val_deg/cos(map.deccen*pi/180);
  
    % now interp sim map to the bin centers of the ra/dec grid as used in
    % reduc_makecomap - the only point of doing this is to strip the
    % area down to the minimum necessary to speed up conv2 and
    % interp2 later
    m=get_map_defn(simopt.type);
    m.T=interp2(map.x_tic,map.y_tic,map.T,m.x_tic',m.y_tic);
    m.Q=interp2(map.x_tic,map.y_tic,map.Q,m.x_tic',m.y_tic);
    m.U=interp2(map.x_tic,map.y_tic,map.U,m.x_tic',m.y_tic);
    % overwrite sim map with re-interpolated version
    map=m;
  
  case 'cel'
    % All sky maps in ra,dec
  
    % load the sim map
    load(simopt.sigmapfilename);
  
    %cel maps might be defined on [-180,180], shift to [0,360]
    if(any(map.x_tic<0))
      res=360/size(map.I,2);
      map.I = circshift(map.I,[0 ceil(180/res)]);
      map.Q = circshift(map.Q,[0 ceil(180/res)]);
      map.U = circshift(map.U,[0 ceil(180/res)]);
      map.x_tic=map.x_tic+180;
    end
  
    % now interp sim map to the bin centers of the ra/dec grid as used in
    % reduc_makecomap 
    m=get_map_defn(simopt.type);
    [xx,yy]=meshgrid(map.x_tic,map.y_tic);
    m.T=interp2(xx,yy,map.I,m.x_tic',m.y_tic);
    m.Q=interp2(xx,yy,map.Q,m.x_tic',m.y_tic);
    m.U=interp2(xx,yy,map.U,m.x_tic',m.y_tic);
    % overwrite sim map with re-interpolated version
    map=m;
  
  case 'comap'
    % when doing leakage sims we want to reload an output map from
    % reduc_makecomap and re-sample from that
    % note that in this case map struct will already be 2 elements long
    load(simopt.sigmapfilename);
    if(simopt.ukpervolt(1)~=1)
      warning('Using comap map and ukpervolt not equal to 1!')
    end
  
  case 'healmap'
    % healpix maps  
    nhmap=length(simopt.sigmapfilename);
    
    % How many support points for interpolating to different beam sizes do we have?
    if ~isempty(simopt.sigmapbeamfwhms)
      nSmooth=length(simopt.sigmapbeamfwhms);
    else
      nSmooth = 1;
    end
    
    % determine the relevant area of sky
    m=get_map_defn(simopt.type);

    % load in a healpix map from B03,WMAP,synfast etc
    for i=1:nhmap
      for k=1:nSmooth
        
        % Generate the filename of the synfast map from the combination
        % of simopt.sigmapfilename and simopt.sigmapbeamfwhms.
        % 
        filename = simopt.sigmapfilename{i};
        if nSmooth>1
          filename = strrep(filename,'.fits','');
          smooth_ext = sprintf('_s%05.2f',simopt.sigmapbeamfwhms(k));
          smooth_ext = strrep(smooth_ext,'.','p');
          l=strfind(filename,'_sxxxx');
          if ~isempty(l)
            % replace _sxxxx string
            le=l+6;
            if le<=numel(filename)
              frest=filename(le:end);
            else
              frest='';
            end
            filename = [filename(1:l-1) smooth_ext frest '.fits'];
          else
            % append to filename
            filename = [filename smooth_ext '.fits'];
          end
        end
          
        % get the healpix map
	disp(sprintf('loading healpix map %s',filename));
        hmap=read_fits_map(filename);
	
        switch hmap.ordering
         case 'RING'
          [theta,phi]=pix2ang_ring(hmap.nside,hmap.pixel);
         case 'NESTED'
          [theta,phi]=pix2ang_nest(hmap.nside,hmap.pixel);
        end
        ra=phi*180/pi; dec=90-theta*180/pi;
        
        % If map is in galactic coords as specified by simopt.coord,
        % rotate ra/dec to galactic as well in order to cut on the right
        % part of the sky:
        if strcmp(simopt.coord,'G')
          [ra,dec]=euler(ra,dec,2,0);
        end
        
        % cut down the map to relevant area of sky
        if(m.hx <180)
          ra(ra>180)=ra(ra>180)-360;
        end
        keepind=find(ra>=m.lx & ra<=m.hx & dec>=m.ly & dec<=m.hy);
        hmap.map=hmap.map(keepind,:);
        hmap.pixel=hmap.pixel(keepind);
        
        % If polarization convention is Healpix change to IAU
        if(strcmp(hmap.polcconv,'COSMO')|strcmp(hmap.polcconv,'UNKNOWN'))
          % flip the sign of U and all its derivatives
          n=3:3:size(hmap.map,2);
          hmap.map(:,n)=-hmap.map(:,n);
          hmap.polcconv='IAU';
          if strcmp(hmap.polcconv,'UNKNOWN')
            warning('Input map has UNKNOWN polcconv - interpreting as healpix convention');
          end
        end

        map(i,k)=hmap;
      end % loop over k
    end % loop over i
   
    if strcmp(simopt.siginterp,'linear') 
      if ~isfield(simopt,'curveskyrescale') || ~simopt.curveskyrescale
        % for simopt.siginterp='linear':
        % make "intermediate map" flat map with specified pixel size
        if size(map,2)>1
          error(['linear interp with beamwidth interpolation '...
                 'is invalid simopt combination']);
        end
        
        for i=1:length(map)
          m.pixsize=simopt.interpix;
          m.nx=round(m.xdos/m.pixsize);
          m.ny=round(m.ydos/m.pixsize);
          sx=(m.hx-m.lx)/m.nx;
          m.x_tic=m.lx+sx/2:sx:m.hx-sx/2;
          sy=(m.hy-m.ly)/m.ny;
          m.y_tic=m.ly+sy/2:sy:m.hy-sy/2;
          
          % go from healpix map to long/lat grid
          if numel(map(1).type)<18
            lmap=healpix_to_longlat(map(i),m.x_tic,m.y_tic,simopt.coord);
          else
            [xx,yy]=meshgrid(m.x_tic,m.y_tic);
            lmap=healpix_interp(map(i),xx(:),yy(:),'taylor',simopt.coord);
            lmap=reshape(lmap,[m.ny,m.nx,3]);
          end
          
          m.T=lmap(:,:,1);
          
          % test if Q/U pol available
          if(numel(map(1).type) > 1) 
            m.Q=lmap(:,:,2);
            m.U=lmap(:,:,3);
          else
            m.Q=zeros(size(m.T));
            m.U=zeros(size(m.T));
          end
          maptmp(i,1)=m;
        end
        
        map=maptmp;
      
        clear hmap
        clear lmap
        
      end
    end
      
    
 otherwise
  
  map=get_map_defn(simopt.type);
  
end

% Modify signal
map=modify_sig(map,simopt.sig);

% inject point sources
switch simopt.ptsrc
  case 'on'
    map=inject_ptsrc(map,simopt);
end
  
if(0)  % plot the maps
  for i=1:length(map)
    figure(i)
    colormap jet
    subplot(3,1,1)
    imagesc(map(i).x_tic,map(i).y_tic,map(i).T); axis xy; ...
    set(gca,'XDir','reverse'); colorbar
    %caxis([-300,300]);
    subplot(3,1,2)
    imagesc(map(i).x_tic,map(i).y_tic,map(i).Q); axis xy; ...
        set(gca,'XDir','reverse'); colorbar
    %caxis([-20,20]);
    subplot(3,1,3)
    imagesc(map(i).x_tic,map(i).y_tic,map(i).U); axis xy; ...
    set(gca,'XDir','reverse'); colorbar
    %caxis([-20,20]);
  end
  drawnow
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map=modify_sig(map,opt)

% BICEP2-style map structure
if isfield(map(1,1),'T')
  for i=1:size(map,1)
    for j=1:size(map,2)
      switch opt

        case 'nopol'
          map(i,j).Q=zeros(size(map(i,j).Q));
          map(i,j).U=zeros(size(map(i,j).U));

        case 'onlypol'
          map(i,j).T=zeros(size(map(i,j).T));

        case 'onlyq'
          map(i,j).T=zeros(size(map(i,j).T));
          map(i,j).U=zeros(size(map(i,j).U));

        case 'test'
          map(i,j).T=zeros(size(map(i,j).T));
          map(i,j).Q=zeros(size(map(i,j).Q));
          map(i,j).U=zeros(size(map(i,j).U));
          % stripe down the middle
          k=abs(map(i,j).x_tic-mean(map(i,j).x_tic))<2;
          map(i,j).T(repmat(k,[size(map(i,j).T,1),1]))=100;
          map(i,j).Q=map(i,j).T;
          map(i,j).U=map(i,j).T;

        case 'none'
          map(i,j).T=zeros(map(i,j).ny,map(i,j).nx);
          map(i,j).Q=zeros(map(i,j).ny,map(i,j).nx);
          map(i,j).U=zeros(map(i,j).ny,map(i,j).nx);

      end
    end
  end

% Healpix data
elseif isfield(map(1,1),'type')
  for i=1:size(map,1)
    for j=1:size(map,2)

      % Identify polarization and Q pol columns by 'type' field
      ispol=false(size(map(i,j).type));
      isq=ispol;
      for k=1:length(ispol)
        ispol(k)=ismember(map(i,j).type{k}(1),'QU') || ismember(map(i,j).type{k}(end),'QU');
        isq(k)=map(i,j).type{k}(1)=='Q' || map(i,j).type{k}(end)=='Q';
      end

      switch opt

        case 'nopol'
          map(i,j).map(:,ispol)=0;

        case 'onlypol'
          map(i,j).map(:,~ispol)=0;

        case 'onlyq'
          map(i,j).map(:,~isq)=0;

        case 'test'
          error('Cannot do taylor interpolation on a test pattern!');

        case 'none'
          map(i,j).map(:,:)=0;

      end
    end
  end

% Some other kind of map?
else
  error(['Unrecognized map type in modify_sig.']);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map=inject_ptsrc(map,simopt)

for j=1:length(simopt.src.name)
  disp(sprintf('Injecting pntsrc %s at RA=%.2f, DEC=%.2f',simopt.src.name{j},simopt.src.ra(j),simopt.src.dec(j))); 
  for i=1:size(map,1)
    for j=1:size(map,2)
      for k=1:size(map,3)
        xi=find(abs(map(i,j,k).x_tic-simopt.src.ra(j))==min(abs(map(i,j,k).x_tic-simopt.src.ra(j))));
        yi=find(abs(map(i,j,k).y_tic-simopt.src.dec(j))==min(abs(map(i,j,k).y_tic-simopt.src.dec(j))));
        map(i,j,k).T(yi,xi)=simopt.src.int(i,j);   
      end
    end
  end
end

return
