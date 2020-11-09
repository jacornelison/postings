function [p,ind]=get_array_info(dat,beamcen,chi,beamwid,diffpoint,flags,epsilon,abscal,polofs)
% [p,ind]=get_array_info(dat,beamcen,chi,beamwid,diffpoint,flags,epsilon,abscal,polofs)
%
% Get array parameters based on date in the tag and modified by the
% various switches
%
% dat = date info required for - e.g. 20120607 numeric or string
% 
% (beamcen, chi, beamwid, diffpoint, epsilon, abscal, polofs) can all be:
%    'ideal'  = values as read from fp_data files
%    filename = values as read from a specific file
%    'obs'    = values as read from designated external files - in
%    this case the requested date will be used to determine which
%    file to read
%    'perazdk' = similar to 'obs' but uses az/dk info as well as the date
%
% "flags" defines additional channel
% flagging in addition to the fp_data file flag column.
% flags is a struture of arrays defining additional datestamped files
% to read and how to cut on the values they contain. e.g: The below
% will cause files aux_data/beams/beams_dummy_rx0_20110101.csv to be
% read and a cut to be made on any channels with fwhm_maj outside the
% range 0.4 to 0.6
%
% flags.filebase={'beams/beams_dummy'}
% flags.par_name={'fwhm_maj'}
% flags.low=.4
% flags.high=.6
%
% [p,ind]=get_array_info(20100106);
% [p,ind]=get_array_info(B3_20160106);

if(~exist('dat','var'))
  dat=[];
end
if(~exist('beamcen','var'))
  beamcen=[];
end
if(~exist('chi','var'))
  chi=[];
end
if(~exist('beamwid','var'))
  beamwid=[];
end
if(~exist('diffpoint','var'))
  diffpoint=[];
end
if(~exist('flags','var'))
  flags=[];
end
if(~exist('epsilon','var'))
  epsilon=[];
end
if(~exist('abscal','var'))
  abscal=[];
end
if(~exist('polofs','var'))
  polofs=[];
end

if(isempty(beamcen))
  beamcen='ideal';
end
if(isempty(chi))
  chi='ideal';
end
if(isempty(beamwid))
  beamwid='ideal';
end
if(isempty(diffpoint))
  diffpoint='ideal';
end
if(isempty(epsilon))
  epsilon='ideal';
end
if(isempty(abscal))
  abscal='ideal';
end
if(isempty(polofs))
  polofs='ideal';
end

% Default date varies by experiment
if isempty(dat)
  expt=get_experiment_name();
  switch(expt)
   case 'bicep3', dat=20160101;
    otherwise, dat=20120106;
  end
end

% expand to allow unified style tag name indicating which expt
expt='';
if(isletter(dat(1)))
  expt=dat(1:2);
  dat=dat(4:end);
end

if(ischar(dat))
  if length(dat)>8
    % should be a tag
    dat=parse_tag(dat);
  else
    dat=str2num(dat);
  end
end

% GET DATA ON EACH FOCAL PLANE FROM THE fp_data FILES

% get the top level file which maps rx to fpu's
map=find_file_by_date(dat,sprintf('aux_data/%s/fp_data/fp_data_master',expt));

% for each listed rx
for i=sort(unique(map.rx))'
  
  % cut down to relevant lines of master file
  mapc=structcut(map,i==map.rx);
  
  % look for line where requested date>start and <end
  after=dat>=mapc.startdate;
  before=dat<mapc.enddate;
  hit=and(after,before);
  if(sum(hit)>1)
    error('Two lines apply to same date in fp_data_master file');
  elseif sum(hit) == 0
    error('Date is outside of all fp_data_master file begin/end dates');
  end

  % get the relevant file for this fpu
  pp=find_file_by_date(dat,sprintf('aux_data/%s/fp_data/fp_data_%s',expt,mapc.fpu{hit}));
  % add a field which specifies which rx this is
  pp.rx=i*ones(size(pp.tile));
  % store angle by which we must rotate this fpu (the drum angle of the cryostat in
  % which it is mounted)
  drumangle(i+1)=mapc.drumangle(hit);
  pp.drumangle=drumangle(i+1)*ones(size(pp.tile));
  
  % allow to to do additional flagging over that done in p.flag
  if(~isempty(flags))
    if ~iscell(flags.filebase)
      flags.filebase={flags.filebase};
      flags.par_name={flags.par_name};
    end
    % To support per-frequency limits, *before* the entire focal plane has
    % been enumerated, allow the bounds to be specified per-receiver. Each
    % row corresponds to a receiver, and replicate if a single row of
    % common thresholds was provided.
    if numel(unique(map.rx)) ~= 1
      if size(flags.low,1)==1
        flags.low = repmat(flags.low, numel(unique(map.rx)), 1);
      end
      if size(flags.high,1)==1
        flags.high = repmat(flags.high, numel(unique(map.rx)), 1);
      end
    end
    for j=1:length(flags.filebase)
      fd=find_file_by_date(dat,sprintf('aux_data/%s/%s_%s',expt,flags.filebase{j},mapc.fpu{hit}));
      fv=getfield(fd,flags.par_name{j});

      % Attempt to use known flags first. Return is uint64, so explicitly cast
      % to double.
      [nextf,err] = conv_fpdata_flags(flags.par_name(j));
      nextf = double(nextf);
      if ~isempty(err)
        warning(sprintf(['Unregistered channel flag ''%s''. Probably ', ...
            'want to add to enum_fpdata_flags()'], flags.par_name{j}))

        knownflags = enum_fpdata_flags();
        maxknown = max(horzcat(knownflags(:).flag));

        % Dynamically assign a bit value greater than any known value,
        % including greater than enumerated list.
        nextf = 2^ceil(log2(max(pp.flag)));
        nextf = max(nextf, double(bitshift(maxknown,1)));
      end
      pp.flag=pp.flag+nextf*(fv<flags.low(i+1,j)|fv>flags.high(i+1,j));
    end
  end

  % create or add to p structure which describes all focal planes
  if(~exist('p','var'))
    p=pp;
  else
    p=structcat(1,[p,pp]);
  end
  
end

% add p.mce=p.rx if mce isn't specified in fp_data
if ~isfield(p,'mce')
  p.mce=p.rx;
end

% add sequential indices for columns and tiles.  These are used in
% various ways in cuts.
p.sequential_mce_col=make_sequential_index(p.rx,p.mce,p.mce_col);
p.sequential_tile=make_sequential_index(p.rx,p.tile);

% put the expt name into p structure
p.expt=get_experiment_name(expt);

% gen lists of channels 
ind=make_ind(p);

% create ideal chi reference theta
p.chi_thetaref=p.theta;

% If not already present get detector parameters from additional file
if(~isfield(p,'r_sh'))
  dp=find_file_by_date(dat,sprintf('aux_data/%s/detparams/detparams',expt));
  for fn=fieldnames(dp)'
    p.(fn{1})=dp.(fn{1});
  end
end

% MODIFY ANY REQUESTED PARAMETERS

% if requested update beam centroid parameters
switch beamcen
  case 'ideal'
    % do nothing - keep what fp_data said
    b=[];
  case 'obs'
    % get what it says in beamcen file
    b=find_file_by_date(dat,sprintf('aux_data/%s/beams/beamcen',expt));
  case 'perazdk'
    if any(strcmp(phase,{'B','E'}))
      az='BE';
    elseif any(strcmp(phase,{'C','F'}))
      az='CF';
    end
    b=find_file_by_date(dat,sprintf('aux_data/%s/beams/beamcen_az%s_dk%03d',expt,az,dk));
  otherwise
    % get a specific file
    b=get_beam_params(beamcen,expt);
end
if ~isempty(b)
  % Ensure no A/B offsets are added here: to make sure,
  % average the A and B beam centers.
  % lat/lon of detectors, possibly with A/B offsets
  [lat,lon]=reckon(0,0,b.r,b.theta+90);
  % distance and angle of offset
  [R,AZ]=distance(lat(ind.a),lon(ind.a),lat(ind.b),lon(ind.b));
  % step back by half of offset
  [lat(ind.a),lon(ind.a)]=reckon(lat(ind.a),lon(ind.a),R/2,AZ);
  lat(ind.b)=lat(ind.a);
  lon(ind.b)=lon(ind.a);

  % go back to r/theta
  p.r=distance(0,0,lat,lon);
  p.theta=azimuth(0,0,lat,lon)-90;

  % Could do same thing using 'paircenter' function
  % p.r=b.r;
  % p.theta=b.theta;
  % [p.r(ind.a),p.theta(ind.a)]=paircenter(b.r(ind.a),b.theta(ind.a),b.r(ind.b),b.theta(ind.b));
  % p.r(ind.b)=p.r(ind.a);
  % p.theta(ind.b)=p.theta(ind.a);
end

% if requested update beamwidths
switch beamwid
  case 'ideal'
    % do nothing - keep what fp_data said
    b=[];    
  case 'obs'
    % get what it says in beamwid file
    b=find_file_by_date(dat,sprintf('aux_data/%s/beams/beamwid',expt));
  case 'zero'
    b.fwhm_maj=zeros(size(p.fwhm_maj));
    b.fwhm_min=b.fwhm_maj;
    b.alpha=b.fwhm_maj;
  otherwise
    % get a specific file
    b=get_beam_params(beamwid,expt);
end
if ~isempty(b)
  if(isfield(b,'fwhm_maj'))
    p.fwhm_maj=b.fwhm_maj;
    p.fwhm_min=b.fwhm_min;
    p.alpha=b.alpha;
  end
  if(isfield(b,'sigma'))
    % convert sigma,c,p to fwhm_maj,fwhm_min,alpha
    % alpha is measured from the respective \vec{e_{theta}} in
    % direction of theta in the same fashion as the polarization
    % angle chi - c and p however are in x'y' coordinates
    [p.fwhm_maj,p.fwhm_min,p.alpha]=egauss2_scp2mmt(b.sigma,b.c,b.p);
    p.alpha=p.alpha-p.theta;
    % in 2014 analysis we added the diff values to the file - not
    % much point - see 20150910_diffellip_refined
    if(isfield(b,'dp'))
      p.dp_fp=b.dp;
      p.dc_fp=b.dc;
    else
      p.dp_fp=nan(size(b.p)); p.dc_fp=p.dp_fp;
      % make sure we don't take the difference of orphan values
      % where one of any pair is zero set the other zero
      c=b;
      c.p(ind.b(c.p(ind.a)==0))=0; c.p(ind.a(c.p(ind.b)==0))=0; 
      c.c(ind.b(c.c(ind.a)==0))=0; c.c(ind.a(c.c(ind.b)==0))=0;
      p.dp_fp(ind.a)=c.p(ind.a)-c.p(ind.b);
      p.dc_fp(ind.a)=c.c(ind.a)-c.c(ind.b);
      p.dp_fp(ind.b)=p.dp_fp(ind.a); p.dc_fp(ind.b)=p.dc_fp(ind.a);
    end
    % We need these in focal plane coords as stored above for the
    % purposes of making diff ellip jacks (when measured in focal
    % plane coords dp means T->E)
    % However for doing diff ellip subtraction we need them in sky
    % coords (at dk=0) so convert to fwhm_maj,fwhm_min,alpha form
    % to make the rotation simple
    avsig=nan(size(b.sigma));
    avsig(ind.a)=(b.sigma(ind.a)+b.sigma(ind.b))/2;
    avsig(ind.b)=avsig(ind.a);
    % convert diff to fwhm_maj,fwhm_min,alpha form
    [p.dfwhm_maj,p.dfwhm_min,p.dalpha]=egauss2_scp2mmt(avsig,p.dc_fp,p.dp_fp);
    p.dalpha=p.dalpha-p.theta;
    % reverse the transform (for testing)
    %[p.s,p.dc,p.dp]=egauss2_mmt2scp(p.dfwhm_maj,p.dfwhm_min,p.dalpha+p.theta);
  end
end

% if requested add in A/B offsets
switch diffpoint
  case 'ideal'
    % do nothing - keep what fp_data said
    b=[];
  case 'obs'
    % get what it says in diffpoint file
    b=find_file_by_date(dat,sprintf('aux_data/%s/beams/diffpoint',expt));
  otherwise
    % get a specific file
    b=get_beam_params(diffpoint,expt);
end
if ~isempty(b)
  % Note this code doesn't give you back exactly the same (r,theta) you put in.

  % lat/lon of detectors with A/B offsets
  [lat,lon]=reckon(0,0,b.r,b.theta+90);
  
  % distance and angle of offset
  [R,AZ]=distance(lat(ind.a),lon(ind.a),lat(ind.b),lon(ind.b));

  % lat/lon of unperturbed detector centroids
  [lat,lon]=reckon(0,0,p.r,p.theta+90);

  % perturb the beam centers to have the same offset magnitude and direction as in the
  % diffpoint beams file
  [lat(ind.a),lon(ind.a)]=reckon(lat(ind.a),lon(ind.a),-R/2,AZ);
  [lat(ind.b),lon(ind.b)]=reckon(lat(ind.b),lon(ind.b),R/2,AZ);
  
  % go back to r/theta
  r=distance(0,0,lat,lon);
  theta=azimuth(0,0,lat,lon)-90;

  % Following version would give you back what you put in:
  %  [rcen,thetacen,xp,yp]=paircenter(b.r(ind.a),b.theta(ind.a),b.r(ind.b),b.theta(ind.b));
  %  r=p.r;
  %  theta=p.theta;
  %  [r(ind.a),theta(ind.a)]=xpyp_to_rtheta(r(ind.a),theta(ind.a),xp{1},yp{1});
  %  [r(ind.b),theta(ind.b)]=xpyp_to_rtheta(r(ind.b),theta(ind.b),xp{2},yp{2});

  isgood=~isnan(r) & ~isnan(theta);
  p.r(isgood)=r(isgood);
  p.theta(isgood)=theta(isgood);
end

% if requested update chi
switch chi
  case 'ideal'
    % do nothing - keep what fp_data said
    qc=[];
  case 'obs'
    % get what it says in chi file
    qc=find_file_by_date(dat,sprintf('aux_data/%s/chi/chi',expt));
  otherwise
    % get a specific file
    qc=ParameterRead(chi);
end
if ~isempty(qc)
  p.chi=qc.chi;
  p.chi_thetaref=qc.theta;
end

% if requested update epsilon
switch epsilon
  case 'ideal'
    % do nothing - keep what fp_data said
    qe=[];
  case 'obs'
    % get what it says in epsilon file
    qe=find_file_by_date(dat,sprintf('aux_data/%s/epsilons/epsilon',expt));
  otherwise
    % get a specific file
    qe=ParameterRead(epsilon);
end
if ~isempty(qe)
  p.epsilon=qe.epsilon;
end

% if requested update abscal
switch abscal
  case 'ideal'
    qa=[];
  case 'obs'
    qa=find_file_by_date(dat,sprintf('aux_data/%s/abscal/abscal',expt));
  otherwise
    qa=ParameterRead(abscal);
end
if ~isempty(qa)
  p.ukpv=qa.ukpv;
else
  % set all channels to nominal overall value
  ukpervolt = get_ukpervolt(num2str(dat),p.expt);
  % expand out to have the number of entries in the p structure
  % assumes that get_ukpervolt delivers the same number of entries
  % per rx - that works for per rx and per tile.
  ukpervolt = repmat(ukpervolt,size(ind.e,2)/numel(ukpervolt),1);
  p.ukpv    = ukpervolt(:);
end

% if requested add in polarization offset ("Roger effect")
switch polofs
  case 'ideal'
    qr=[];
  case 'obs'
    qr=find_file_by_date(dat,sprintf('aux_data/%s/polofs/polofs',expt));
  otherwise
    qr=ParameterRead(polofs);
end
if ~isempty(qr)
  p.polofs_x=qr.polofs_x;
  p.polofs_y=qr.polofs_y;
end

% rotate quantities by the drum angle
for i=1:numel(drumangle)
  rotind=p.rx==i-1;
  p.theta(rotind)=p.theta(rotind)-drumangle(i);
  p.chi_thetaref(rotind)=p.chi_thetaref(rotind)-drumangle(i);
  if isfield(p,'polofs_x')
    tmp=p.polofs_x(rotind)*cosd(drumangle(i))-p.polofs_y(rotind)*sind(drumangle(i));
    p.polofs_y(rotind)=p.polofs_y(rotind)*cosd(drumangle(i))+p.polofs_x(rotind)*sind(drumangle(i));
    p.polofs_x(rotind)=tmp;
  end
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b=get_beam_params(type,expt)
% read beam parameters from file allowing latitude in name
% specification

[fdir fname fext]=fileparts(type);
if isempty(fdir)
  fdir=sprintf('aux_data/%s/beams',expt);
end
if isempty(fext)
  fext='.csv';
end

b=ParameterRead(fullfile(fdir,[fname fext]));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function i=make_sequential_index(varargin)
key=[];
for j=1:length(varargin)
  key=[key,reshape(varargin{j},[],1)];
end
% Inf replaces NaN because every NaN is a special snowflake
badval=any(~isfinite(key),2);
key(badval,:)=Inf;
[tmp j i]=unique(key,'rows');
i(badval)=NaN;

return

