function bm_timestream(bmopt)
%
% function bm_timestream(bmopt)
% 
% TSG 2018-06-27
%
% First step in processing raw beam map data.
% For a given FFBM schedule, load up arcfiles one at a time and:
%  -- Remove/append incomplete halfscans, so each file has an 
%     integer number of halfscans
%  -- get pointing model and invert -> hor coords
%  -- deconvolve detector timestreams with pipeline standard kernel
%  -- swap in for simulated data if desired 
%  -- demodulated and downsample
%  -- get coordinate timestreams in xp/yp, x/y, and apparent az/el
%  -- apply masking at timestream level (ground/mast/SPT)
%  -- save output, in ~arcfile-size units
% 
% INPUTS (should all be passed in with bmopt)
%
%   expt:     'bicep2','keck','bicep3'
%   year:     Bicep2: 2012
%             Keck: 2012+
%             Bicep3: 2016+
%   number:   which run # from this FFBM campaign? From runfile
%   run:      BICEP2 runs: 'r5','r8'
%             Keck runs:   'k0r4', 'k0r4t2' (tile 2 only), highbay
%                          'k083' = k7 run 2013-10, highbay
%                          'k201?' = Keck ffbm @ Pole (12,13,14,15)
%                          'k201?fsl' = Keck sidelobe @ Pole (13,14)
%             BICEP3 runs: 'hb3r5' = B3 run 5, highbay
%                          'b3r6' = B3 run 6 @ Pole, 2015-02
%                          'b3r8' = B3 run 8 @ Pole, 2016-02
%                          'b3r9' = B3 run 9 @ Pole, 2017/18-02
%   source_az:     source type CENTER azimuth position (deg)
%   source_el:     source type CENTER el position (deg)
%
% OPTIONAL
%
%  notlastscan:  = 0 (default) demodulates one file at a time
%                = 1 (recommended when analyzing entire schedule)
%                    keeps the last scan if it's incomplete, 
%                    does not demodulate that scan, and
%                    attaches it to the following file.  Must start from
%                    the beginning of the beam map schedule when
%                    using this option.  Slows computation time but
%                    ensures no data in the middle of a scan is 
%                    "thrown out".
%  demodtype:     'square' (default), 'choplet'
%                  Choplet does not work as of this pipline
%                  upgrade, need to make it compatible with
%                  cleaned chop ref.
%  cleanref:     = 0 demodulates using chop ref as-is from PMAC
%                = 1 (default) puts chop ref through cleanchopref.m and runs
%                    square_demod with highres option.
%                    Recommended due to pathologies in most of the 
%                    real chop refs from FFBM.
%  demodshift:   number of samples by which to shift ref signal - we have
%                hard-coded default values for various times, but can 
%                change this if we want.  Gets saved to demodded
%                file.
%  margin:       For timestream mirror masking:
%                We don't want to mask at exactly the mirror dimensions
%                since we know there are beam distortions near edge of 
%                mirror.  This number (deg) is how much to cut into the
%                mask at each edge.  Default = 2.5 deg.
%  btsimopt:     = 0 to keep real data (default)
%                If nonzero, should be struct containing information about
%                BM timestream sim to use:
%           .mapopt: filename string of input map, or struct describing
%                    arbitrary map (see bmtodsim_inputmap)
%           .refopt: describes how to modulate simulated signal 
%                    (see bmtodsim_insertsimdata)
%  pairdiffsum:    = 0 do individual detector timestreams (default)
%                  = 1 make pair sum and pair diff timestreams,
%                    taking (A+B)/2 and (A-B)/2 pre-demodding.
%                    In this case, need fit amplitudes from
%                    individual detector beams.
%  beamparammatdir:  directory to beam param .mat file (pairdiffsum=1)
%                    default 'beammaps/beamparams/'
%  beamparammatfile: name of beam param .mat file (pairdiffsum=1)
%                    default 'beamparams_year.mat'
%  timestreamdir:  directory in which to save output timestream
%                  default is 'beammaps/timestream/'
%  update:       option to not reprocess already existing data products
%                 = 0 (default) will re-process all files
%                 = 1 will only process missing files


% Parse input and set defaults.  Use get_bm_info.
expt = bmopt.expt;
year = bmopt.year;
number = bmopt.number;
run = bmopt.run;
source_az = bmopt.source_az;
source_el = bmopt.source_el;

pp = get_bm_info(year);
t1 = pp.t1{number};
t2 = pp.t2{number};
sched = pp.sched{number};
if strcmp(expt,'keck')
  mirror = pp.mirror{number};
else
  mirror = '';
end
if ~isfield(bmopt,'notlastscan')
  notlastscan = 0;
  bmopt.notlastscan = notlastscan;
else
  notlastscan = bmopt.notlastscan;
end
if ~isfield(bmopt,'demodtype');
  demodtype = 'square';
  bmopt.demodtype = demodtype;
else
  demodtype = bmopt.demodtype;
end
if ~isfield(bmopt,'timestreamdir')
  timestreamdir = 'beammaps/timestream';
  bmopt.timestreamdir = timestreamdir;
else
  timestreamdir = bmopt.timestreamdir;
end
if ~isfield(bmopt,'cleanref')
  cleanref = 1;
  bmopt.cleanref = cleanref;
else
  cleanref = bmopt.cleanref;
end
if ~isfield(bmopt,'demodshift')
  demodshift = [];
else
  demodshift = bmopt.demodshift;
end
if ~isfield(bmopt,'margin')
  margin = 2.5;
else
  margin = bmopt.margin;
end
if ~isfield(bmopt,'btsimopt')
  btsimopt = 0;
  bmopt.btsimopt = btsimopt;
else
  btsimopt = bmopt.btsimopt;
end
if ~isfield(bmopt,'pairdiffsum')
  pairdiffsum = 0;
  bmopt.pairdiffsum = pairdiffsum;
else
  pairdiffsum = bmopt.pairdiffsum;
end
if ~isfield(bmopt,'beamparammatdir')
  beamparammatdir = 'beammaps/beamparams/';
  bmopt.beamparammatdir = beamparammatdir;
else
  beamparammatdir = bmopt.beamparammatdir;
end
if ~isfield(bmopt,'beamparammatfile')
  beamparammatfile = ['beamparams_' num2str(year)];
  bmopt.beamparammatfile = beamparammatfile;
else
  beamparammatfile = bmopt.beamparammatfile;
end
if ~isfield(bmopt,'update')
  update = 0;
  bmopt.update = update;
else
  update = bmopt.update;
end

% Get a list of the relevant files 
files = list_arc_files('arc/',t1,t2);

[p ind] = get_array_info(files{1}(1:8),'ideal');

% For concatenating previous files
e = [];

bmopt
tic;

% If we're using sim data, load the input map now
if isstruct(btsimopt)
  disp('Loading input map for simulated timestreams...')
  [xp_sim yp_sim mapsim] = bmtodsim_inputmap(btsimopt.mapopt);
  btsimopt.xp_sim = xp_sim;
  btsimopt.yp_sim = yp_sim;
  btsimopt.mapsim = mapsim;
  disp(['Done loading input map (' num2str(toc) ' sec).'])
end

% If we're pair summing/differencing, get fit amplitudes now
if pairdiffsum
    w = load([beamparammatdir '/' beamparammatfile]);
    amp = w.beamparam.amp;
end

% Loop over arcfiles
for fl = 1:length(files)
  
  % if update, check if timestream file already exists, don't overwrite
  if update == 1
    flstr = [timestreamdir '/bm_' files{fl}(1,1:15) '.mat']; 
    if exist_file(flstr) 
      disp(sprintf('file %s already exists, skipping to next',flstr));
      continue
    end
  end
    
  % Load the data
  % In rare occasion where the file is corrupted, clear any
  % data being carried over from notlastscan and continue
  % to the next file.
  load_file = char(strcat('arc/',files(fl)));
  disp(['Loading file: ' load_file]);
  try
    d = load_arc(load_file);
    disp(['File loaded (' num2str(toc) ' sec)']);
  catch
    disp(['WARNING: file did not load.  Check if corrupted. ' ...
          'Moving to next arcfile.'])
    e = [];
    continue
  end
  
  % Notlastscan -- if file ends during a halfscan, keep that
  % scan, don't demodulate it, append it to beginning of next file.
  % We do this because square_demod returns junk/NaN at very
  % beginning and end of whatever segment we are demodulating,
  % so we try to make those happen at the turnarounds.
  if notlastscan
    dtmp = d;
    if ~isempty(e)
      clear d;    
      d = structcat(1,[e,dtmp]);
    end
    clear dtmp;
    dtmp = d;
  end
  
  % Check to see that there is valid data
  mark = d.array.frame.features;
  dm = d.mce0.rc1.data_mode(logical(mark));
  if numel(dm) == 0
    disp(['File ' load_file ' is empty'])
    continue 
  end

  % Skip this file if there are no field scans
  fs = get_fieldscans(d);
  if isempty(fs.s)
    disp(['No field scans found in this file.  Moving to next ' ...
          'file.'])
    continue
  end
  
  % Concatenate the 'mceX' fields together
  d = concatmce(d);

  % Deconvolve
  d = get_xferfunc(d,fs,0,0); 
  d = deconv_scans(d,p,1);
  disp(['Deconv scans done (' num2str(toc) ' seconds)']);
  
  % Generate pointing model
  mjd = date2mjd(str2num(files{1}(1:4)),str2num(files{1}(5:6)),...
                 str2num(files{1}(7:8)));
  bmopt.pm = get_pointing_model(mjd,0,d,fs);
  bmopt.pm.az_tilt_ha = 0;
  bmopt.pm.az_tilt_lat = 0;
  bmopt.pm.el_tilt = 0;
  
  % Invert pm to hor coords
  d = get_horcoords(d,bmopt.pm,fs,run);

  % Trim off fields not needed for later analysis
  d = rename_fields(d,p,expt,run);
  
  if notlastscan
    if ~isempty(e)
      % If there is a file before, take the last scan from the file before
      tmp1 = max(find(fs.sf < length(e.mce0.data.fb(:,1))));
    else
      % Else start from the first scan
      tmp1 = 1;
    end
    if fs.ef(end) == length(d.fb(:,1));
      % If the scan goes to the end of the file, 
      % it means that the last scan ends in the 
      % file that comes after this, so skip the last scan
      % save the last bit of the file to go onto the next file
      tmp2 = length(fs.ef) - 1;
      clear e;
      cutscan = structcut(fs,tmp2+1);
      e = cutdstruct(dtmp,cutscan);
      clear dtmp;
    else
      % If the scan doesn't go to the end of the file,
      % don't save the file -- demodulate to end of file
      clear e;
      e = [];
      tmp2 = length(fs.ef);
      clear dtmp;
    end
    if tmp2 > tmp1 & tmp2 ~= 0
      % Keep the sample before the start of the last az scan
      indexcut = [fs.sf(tmp1):fs.ef(tmp2)];
      % Downsampling is messy when we have a struct in a struct
      pointing = d.pointing;
      d = rmfield(d,'pointing');
      d = structcut(d,indexcut);
      pointing.hor = structcut(pointing.hor,indexcut);
      d.pointing = pointing;
      % Also have to cut fs
      fstmp = fs;
      fs = structcut(fstmp,[tmp1:tmp2]);
      fs.sf = fs.sf - fstmp.sf(tmp1) + 1;
      fs.ef = fs.ef - fstmp.sf(tmp1) + 1;
    else 
      continue;
    end
    disp(['Notlastscan done (' num2str(toc) ' seconds)']);
  end % notlastscan
  
  % Swap in simulated timestream data if requested
  if isstruct(btsimopt)
    disp(['Swapping in simulated timestreams...'])
    btsimopt.source_az = source_az;
    btsimopt.source_el = source_el;
    btsimopt.run = run;
    btsimopt.p = p;
    d = bmtodsim_insertsimdata(d,btsimopt);
    demodshift = 0;
    disp(['Sim data done (' num2str(toc) ' sec).'])
  end
  
  % Make chop ref binary
  switch expt
    case 'bicep2'
      sqw = d.ref < mean(d.ref);
    case {'keck','bicep3'}
      sqw = d.ref > mean(d.ref);     
  end
  
  % If B2r8 or K0r4, symmetrically high-pass filter the fb:
  switch run
    case {'r8','k0r4'}
      % Define high-pass filter coefficients:
      [bf af] = butter(1,1/4000,'high');
      d.fb = filtfilt(bf,af,d.fb);
  end

  % Deglitch chop reference for bicep2 data @ pole, summer 2010-11
  switch run
    case {'r81','r82','r83','r84','r85','r86'}
      sqw = deglitch_sqw(sqw);
  end
  
  % Make pair-diff and pair-sum timestreams if requested
  if pairdiffsum
    % Use apprpriate normalization
    A_amp = repmat(amp(ind.la,number)',[size(d.fb,1),1]);
    B_amp = repmat(amp(ind.lb,number)',[size(d.fb,1),1]);
    A = d.fb(:,ind.la) ./ A_amp;
    B = d.fb(:,ind.lb) ./ B_amp;
    d.fb_sum = (A+B) / 2;
    d.fb_diff = (A-B) / 2;
    d = rmfield(d,'fb');
  end
  
  % Demodulate!
  % On rare occasions, this step fails due to large anomalies in 
  % the chop ref.  If that happens, clear data being carried over
  % from notlast scan and move to next file.
  disp('Demodulation start');
  %  try
    d = demod_blocks(d,fs,sqw,run,demodtype,expt,cleanref,...
                     pairdiffsum,demodshift);
    disp(['Demodulation done (' num2str(toc) ' seconds)']);
    %  catch
    %    disp(['WARNING: demod step failed. Check the chop reference. ' ...
    %          'Moving to next arcfile.'])
    %    e = [];
    %    continue
    %  end
  
  % Get the coordinate timestreams in apparent az/el, x/y, xp/yp
  % Have to do this after downsampling, otherwise too
  % computationally expensive.
  disp(['Calculating coordinate systems...']);
  d = get_bmcoords(d,p,ind,bmopt.pm,fs,source_az,source_el,run,pairdiffsum);
  disp(['Coords done (' num2str(toc) ' seconds)']);
  
  % Masking
  disp(['Creating masks...']);
  % Want a (N_index x N_det) matrix mask (1 good, 0 bad) each
  % for ground+SPT+mask 
  d.maskground = get_groundmask(d,p,ind,expt,source_az,source_el);
  % Mask mirror is evaluated using x/y coords, so mask is 
  % (N_index x N_rx).  Returns all 1s for non-Keck.
  d.maskmirror = get_mirrormask(d,p,expt,mirror,sched,margin);
  disp(['Masking done (' num2str(toc) ' seconds)']); 
  
  % Save!
  % If there is a user-input demod shift (not from sims), 
  % save it in the filename
  if ~isempty(demodshift) && ~isstruct(btsimopt)
    demodstr = ['_shift' int2str(demodshift)];
  else
    demodstr = [];
  end
  if pairdiffsum
    pairstr = 'pair_';
  else
    pairstr = '';
  end
  flstr = ...
      [timestreamdir '/bm_' pairstr files{fl}(1,1:15) '_' demodtype demodstr]; 
  disp(['Saving ' flstr])
  save(flstr,'d','-v7.3')
  
end % arcfile loop

disp(['Exiting bm_timestream']);

return






%%%%%%%%%%%%%%%%%%%%%% SUB-FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%

function d = concatmce(d)

% Look for any fields 'mceX' where X > 0 and concat them to mce0
fnames = fieldnames(d);
if length(fnames) > 3
  for ii = 4:length(fnames)
    m = fnames{ii};
    d.mce0 = structcat(2,[d.mce0,getfield(d,m)]);
    d = rmfield(d,m);
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fs = get_fieldscans(d)
  
% Generate field scans struct
% Fast/slow ratio
sampratio = ...
    length(d.antenna0.pmac.fast_az_pos)/length(d.antenna0.pmac.az_pos); 
try
  % When is feature bit 1 on?
  fs = find_scans(d.antenna0.tracker.scan_off(:,1),...
                  bitand(d.array.frame.features,2^1),sampratio,'keepend');
catch
  fs.s = [];
  fs.e = [];
end
disp(['Found scans (' num2str(toc) ' seconds)']);
% Be sure there isn't a crazy turnaround right at the beginning where the
% start of the scan happens after the beginning:
fs = structcut(fs,fs.s < fs.e);
% Also check that there are enough samples between the start/end of scan
% Default Keck xfer function is 22 samples long, give a little buffer
% B3 xfer function may not be correct as is, check this step in future
fs = structcut(fs,fs.e - fs.s > 30);  
for ii = 1:length(fs.e)
  if fs.e(ii) > length(d.mce0.cc.data_rate) 
    fs.e(ii) = length(d.mce0.cc.data_rate); 
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = get_horcoords(d,pm,fs,run)

% First invert pointing model to get horizon coords
switch run
  case {'k2012','k2013','k2014','k2015','k2016','k2017','k2018'} 
    d = invpointing_model(d,pm,0,fs);
  case {'k2013fsl','k2014fsl'}   % Sidelobe mapping
    d = invpointing_model(d,pm);
  case {'b3r6','b3r8','b3r9'} 
    % keck_beam_map_pointing wants where telescope was truly pointing
    % during beam mapping, so we have to undo the changes made
    % to the encoder during the FFBM schedule (only done for B3).
    % Set all encoder_sign to +1 and flip sign of dk encoder offset.
    % Then DON'T use mirror flag on invpointing_model, makes sure
    % no extra sign flip is applied.
    d.antenna0.tracker.encoder_sign = abs(d.antenna0.tracker.encoder_sign);
    d.antenna0.tracker.encoder_off(:,3) = ...
        -d.antenna0.tracker.encoder_off(:,3);
    d = invpointing_model(d,pm,0,fs);
  otherwise                      
    d = invpointing_model(d,pm,1,fs);
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = cutdstruct(d,cutscan)

e = [];

thing1 = fieldnames(d);

for ii = 1:length(thing1)
  thing2 = fieldnames(getfield(d,thing1{ii}));
  for jj = 1:length(thing2)
    thing3 = fieldnames(getfield(d,thing1{ii},thing2{jj}));
    for kk = 1:length(thing3)
      tmp = getfield(d,thing1{ii},thing2{jj},thing3{kk});
      if size(tmp,1) == 0
        e = setfield(e,thing1{ii},thing2{jj},thing3{kk},tmp);
      elseif size(tmp,1) < 4001
        tmpnew = tmp(cutscan.s:cutscan.e,:);
        e = setfield(e,thing1{ii},thing2{jj},thing3{kk},tmpnew);
      elseif size(tmp,1) > 4000
        tmpnew = tmp(cutscan.sf:cutscan.ef,:);
        e = setfield(e,thing1{ii},thing2{jj},thing3{kk},tmpnew);
      else
      end
      clear tmp;
      clear tmpnew;
    end
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = demod_blocks(d,fs,sqw,run,demodtype,expt,cleanref,...
                          pairdiffsum,demodshift)

switch demodtype
  case 'choplet'
    %d.mce0.data.fb = d.fb;
    % What circshift do we want?
    % 2012-02
    %if d.t(1) < date2mjd(2013) 
    %  shiftme = 2;
    % 2013-02
    %elseif d.t(1) < date2mjd(2013,2,14) 
    %  shiftme = 7;
    % 2013-02 to 2013-03
    %elseif d.t(1) > date2mjd(2013,2,18) & d.t(1) < date2mjd(2013,3,10) 
    %  shiftme = 2;
    %end
    % Apply the shift
    %if ~isempty(demodshift)
    %  shiftme = demodshift;
    %end
    %d.antenna0.pmac.fast_aux_input(:,1) = circshift(d.ref,[shiftme,0]);    
    %d.t = utc2datenum(cat(2,d.t,d.t2));
    %d = demod_choplet(d);
    %d = rmfield(d,{'antenna0','fb','mce0','t','t2'})  
    error('demod_choplet not optimized for cleaned chop ref yet')
    
  case 'square'
  
    d.t = d.t(1);

    % Punting on B2 demod shifting...way too complicated
    switch expt
      case 'bicep2'
	shiftme = 0;
	if ~isempty(demodshift)
	  shiftme = demodshift;
	end
      case 'keck'
	switch run
	  case 'h083'
	    shiftme = -1;
	  case {'k2012','k2013','k2013fsl','k2014','k2014fsl','k2015','k2016','k2017','k2018'}
	    if d.t < date2mjd(2013) % 2012
	      shiftme = 1;          % COS max(used to be 2)
	    elseif d.t < date2mjd(2013,2,14) % 2013-02 round 1
	      shiftme = 6;                   % COS max (used to be 7)
	    elseif d.t > date2mjd(2013,2,18) & ...
		  d.t < date2mjd(2013,3,10)  % 2013-02 to 2013-03
	      shiftme = 2;                   % COS max
	    elseif d.t > date2mjd(2014,1,30) & ...
		  d.t < date2mjd(2014,2,14,02,16,04) % 2014-02a
	      shiftme = 4;
	    elseif d.t > date2mjd(2014,2,14,02,16,04) & ...
		  d.t < date2mjd(2014,3,1,07,35,37)  % 2014-02b
	      shiftme = 3;                           % COS max
	    elseif d.t > date2mjd(2014,3,1,07,35,37) & ...
		  d.t < date2mjd(2014,3,13) % 2014-03
	      shiftme = -5;                 % COS max
	    elseif d.t < date2mjd(2015) % 2014-01 FSL mapping on MAPO
	      shiftme = -3;
	    elseif d.t < date2mjd(2016) % 2015-02
	      shiftme = 5;              % COS max (used to be 3)
	    elseif d.t < date2mjd(2016,2,27) % 2016-02 part 1
	      shiftme = 0;                   % COS max (used to be -2)
            elseif d.t < date2mjd(2017) % 2016-02 to 2016-03
              shiftme = -1;             % COS max
            elseif d.t < date2mjd(2018) % 2017-02
              shiftme = 7;              % COS max (used to be 6)
	    elseif d.t < date2mjd(2018,2,7) % 2018 BSNS as source
              %shiftme = 10;                % SIN max
	      shiftme = 3;                  % COS max
	    elseif d.t < date2mjd(2019) % 2018 chopper as source
	      %shiftme = 3;             % SIN max
	      shiftme = 7;              % COS max
	    else
	      disp('No demod phase chosen for this date!  Defaulting to 0');
	      shiftme = 0;
	    end
	end
	if ~isempty(demodshift)
	  shiftme = demodshift;
	end
      case 'bicep3'
	switch run
	  case 'hb3r5'
	    shiftme = -2;
	  case 'b3r6' 
	    if d.t < date2mjd(2015,3,16,23,45,00) % 2015-02 24in chopper
	      shiftme = -3;                       % COS max
	    elseif d.t > date2mjd(2015,3,16,23,45,00) & ...
		  d.t < date2mjd(2016)            % 2015-02 18in chopper
	      shiftme = 3;                        % COS max
	    end
	  case 'b3r8'     % 2016-02
	    %shiftme = 2;  % COS max (before adding deconv)
	    shiftme = 1;  % COS max
	  case 'b3r9'
	    if d.t < date2mjd(2018)   % 2017-02
	      %shiftme = -1;          % COS max (before adding deconv)
	      shiftme = -2;           % COS max
	    elseif d.t < date2mjd(2018,2,5) % 2018 BSNS as source
	      %shiftme = -1;                % SIN max
	      %shiftme = -8;                % COS max (before adding deconv)
	      shiftme = -9;                 % COS max
	    elseif d.t < date2mjd(2019) % 2018 chopper as source
	      %shiftme = -3;            % SIN max
	      %shiftme = -1;            % COS max (before adding deconv)
	      shiftme = -2;             % COS max
	    end
	end
	if ~isempty(demodshift)
	  shiftme = demodshift;
	end
    end
    
    d = rmfield(d,{'t','t2'});

    disp(['Shifting chop ref by ' num2str(shiftme) ' indices.'])
    d.sqw = shiftchopref(sqw,shiftme);  
    
    % Demodulate either per-det timestreams or pairsum/pairdiffs
    % Clean chop ref if requested, then demodulate, then clean
    % up indices and downsample
    if pairdiffsum
      if cleanref == 1
        disp(['Cleaning chop ref...'])
        ch = cleanchopref(d.sqw);
        [e.cos_sum,e.sin_sum,e.ind] = ...
            square_demod(d.fb_sum,ch.refclean_hires,'highres');
        [e.cos_diff,e.sin_diff,e.ind] = ...
            square_demod(d.fb_diff,ch.refclean_hires,'highres');
      else
        [e.cos_sum,e.sin_sum,e.ind] = ...
            square_demod(d.fb_sum,d.sqw);
        [e.cos_diff,e.sin_diff,e.ind] = ...
            square_demod(d.fb_diff,d.sqw);
      end

      % Deal with case where sin/cos components have different length  
      if length(e.sin_sum(:,1)) > length(e.cos_sum(:,1))
        e.sin_sum = e.sin_sum(2:length(e.cos_sum(:,1)) + 1,:);
      elseif length(e.sin_sum(:,1)) < length(e.cos_sum(:,1))
        e.sin_sum = ...
            [e.sin_sum; zeros(length(e.cos_sum(:,1)) - ...
                              length(e.sin_sum(:,1)),size(e.sin_sum,2))];
      end
      if length(e.sin_diff(:,1)) > length(e.cos_diff(:,1))
        e.sin_diff = e.sin_diff(2:length(e.cos_diff(:,1)) + 1,:);
      elseif length(e.sin_diff(:,1)) < length(e.cos_diff(:,1))
        e.sin_diff = ...
            [e.sin_diff; zeros(length(e.cos_diff(:,1)) - ...
                              length(e.sin_diff(:,1)),size(e.sin_diff,2))];
      end
      
      % Downselect indices that are within the field scans
      cc = [];
      for ii = 1:length(fs.sf)
        keep = e.ind(e.ind > fs.sf(ii) & e.ind < fs.ef(ii));
        cc = [cc; keep];
      end
      switch expt
        case {'keck','bicep3'}
          e.cos_sum = e.cos_sum(ismember(e.ind,cc),:); 
          e.sin_sum = e.sin_sum(ismember(e.ind,cc),:); 
          e.cos_diff = e.cos_diff(ismember(e.ind,cc),:); 
          e.sin_diff = e.sin_diff(ismember(e.ind,cc),:);
          e.ind = e.ind(ismember(e.ind,cc),:); 
        case 'bicep2'
          if size(e.cos_sum) == size(e.sin_sum)
            e = structcut(e,ismember(e.ind,cc));
          else
            e.sin_sum = e.sin_sum(1:length(e.cos_sum),:);
            e.cos_sum = e.cos_sum(ismember(e.ind,cc),:);
            e.sin_sum = e.sin_sum(ismember(e.ind,cc),:); 
            e.ind = e.ind(ismember(e.ind,cc),:); 
          end
          if size(e.cos_diff) == size(e.sin_diff)
            e = structcut(e,ismember(e.ind,cc));
          else
            e.sin_diff = e.sin_diff(1:length(e.cos_diff),:);
            e.cos_diff = e.cos_diff(ismember(e.ind,cc),:);
            e.sin_diff = e.sin_diff(ismember(e.ind,cc),:); 
            e.ind = e.ind(ismember(e.ind,cc),:); 
          end

      end
      
      % Downsample the d struct
      pointing = d.pointing;
      d = rmfield(d,'pointing');
      d = structcut(d,cc);
      pointing.hor = structcut(pointing.hor,cc);
      % And put pointing back
      d.pointing = pointing;
      
      d.cos_sum = e.cos_sum; 
      d.sin_sum = e.sin_sum; 
      d.cos_diff = e.cos_diff; 
      d.sin_diff = e.sin_diff; 
      d.ind = e.ind;
      d.demodshift = shiftme;
      d = rmfield(d,'sqw');

    else % per det, not pairsum/diff
      
      if cleanref == 1
        disp(['Cleaning chop ref...'])
        ch = cleanchopref(d.sqw);
        [e.cos,e.sin,e.ind] = square_demod(d.fb,ch.refclean_hires,'highres');
      else
        [e.cos,e.sin,e.ind] = square_demod(d.fb,d.sqw);
      end
    
      % Deal with case where sin/cos components have different length  
      if length(e.sin(:,1)) > length(e.cos(:,1))
        e.sin = e.sin(2:length(e.cos(:,1)) + 1,:);
      elseif length(e.sin(:,1)) < length(e.cos(:,1))
        e.sin = [e.sin; ...
                 zeros(length(e.cos(:,1)) - length(e.sin(:,1)),size(e.sin,2))];
      end
       
      % Downselect indices that are within the field scans
      cc = [];
      for ii = 1:length(fs.sf)
        keep = e.ind(e.ind > fs.sf(ii) & e.ind < fs.ef(ii));
        cc = [cc; keep];
      end
      switch expt
        case {'keck','bicep3'}
          e.cos = e.cos(ismember(e.ind,cc),:); 
          e.sin = e.sin(ismember(e.ind,cc),:); 
          e.ind = e.ind(ismember(e.ind,cc),:); 
        case {'bicep2'}
          if size(e.cos) == size(e.sin)
            e = structcut(e,ismember(e.ind,cc));
          else
            e.sin = e.sin(1:length(e.cos),:);
            e.cos = e.cos(ismember(e.ind,cc),:); %for bicep2 summer 2012 
            e.sin = e.sin(ismember(e.ind,cc),:); 
            e.ind = e.ind(ismember(e.ind,cc),:); 
          end
      end
      
      % Downsample the d struct
      pointing = d.pointing;
      d = rmfield(d,'pointing');
      d = structcut(d,cc);
      pointing.hor = structcut(pointing.hor,cc);
      % And put pointing back
      d.pointing = pointing;
      
      d.cos = e.cos; 
      d.sin = e.sin; 
      d.ind = e.ind;
      d.demodshift = shiftme;
      d = rmfield(d,'sqw');
      
    end

end % Standard demod

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = get_bmcoords(d,p,ind,pm,fs,source_az,source_el,run,pairdiffsum)

% Convert hor coords from double->single to trim file size
% Will do this for all other coords too
d.pointing.az = single(d.pointing.hor.az);
d.pointing.el = single(d.pointing.hor.el);
d.pointing.dk = single(d.pointing.hor.dk);
d.pointing = rmfield(d.pointing,'hor');

% Need mirror/mount parameters for other coord systems
[mount mirror source] = get_mms_args(run);
source.azimuth = source_az; % degrees
source.height = tan(source_el * pi/180) * source.distance; % meters

% Apparent az/el 
if pairdiffsum
  ptemp.r = p.r(ind.la);
  ptemp.theta = p.theta(ind.la);
  ptemp.drumangle = p.drumangle(ind.la);
  ptemp.chi = p.chi(ind.la);
  ptemp.chi_thetaref = p.chi_thetaref(ind.la);
else
  ptemp = p;
end
[az_ap, el_ap, pa_out] = keck_beam_map_pointing(d.pointing.az,...
                           d.pointing.el,d.pointing.dk,mount,...
                           mirror,source,ptemp,'Legacy');
az_ap(az_ap > 90) = az_ap(az_ap > 90) - 360;
az_ap(az_ap < -270) = az_ap(az_ap < -270) + 360;
d.pointing.az_ap = single(az_ap);
d.pointing.el_ap = single(el_ap);

% x/y
% Only need one timestream per rx
% Plotting: for  both x/y and xp/yp, when plotting you want to do:
%  >> imagescnan(y,x,map);
%  >> set(gca,'YDir','normal','XDir','reverse')
% to get standard parity maps.
ptemp.drumangle = flipud(unique(p.drumangle));
ptemp.r = zeros(size(ptemp.drumangle));
ptemp.theta = zeros(size(ptemp.drumangle));
ptemp.chi = zeros(size(ptemp.drumangle));
ptemp.chi_thetaref = zeros(size(ptemp.drumangle));
[r,theta,psi] = keck_beam_map_pointing(d.pointing.az,...
    d.pointing.el,d.pointing.dk,mount,mirror,source,ptemp);
d.pointing.y = ...
    single(r .* sin((theta+repmat(ptemp.drumangle,1,size(r,1))')* pi/180));
d.pointing.x = ...
    single(r .* cos((theta+repmat(ptemp.drumangle,1,size(r,1))')* pi/180));

% xp/yp
% Due to distortion effects we shouldn't just translate x/y to xp/yp.
if pairdiffsum
  ptemp.r = p.r(ind.la);
  ptemp.theta = p.theta(ind.la);
  ptemp.drumangle = p.drumangle(ind.la);
  ptemp.chi = p.chi(ind.la);
  ptemp.chi_thetaref = p.chi_thetaref(ind.la);
else
  ptemp = p;
end
[r,theta,psi] = keck_beam_map_pointing(d.pointing.az,...
    d.pointing.el,d.pointing.dk,mount,mirror,source,ptemp);
d.pointing.yp = ...
    single(r .* sin((theta+repmat(ptemp.drumangle,1,size(r,1))')* pi/180));
d.pointing.xp = ...
    single(r .* cos((theta+repmat(ptemp.drumangle,1,size(r,1))')* pi/180));

% xp/yp timestream for buddy map
% Again, can't just translate from standard xp/yp due to distortion
if pairdiffsum
  pbuddy.r = p.r(ind.la);
  pbuddy.theta = p.theta(ind.la);
  pbuddy.drumangle = p.drumangle(ind.la);
  pbuddy.chi = p.chi(ind.la);
  pbuddy.chi_thetaref = p.chi_thetaref(ind.la);
else
  pbuddy = p;
  pbuddy.theta = pbuddy.theta - 180;
end
pbuddy = p;
pbuddy.theta = pbuddy.theta - 180;
[rb,thetab,psib] = keck_beam_map_pointing(d.pointing.az,...
    d.pointing.el,d.pointing.dk,mount,mirror,source,pbuddy);
d.pointing.yp_buddy = ...
    single(rb .* sin((thetab+repmat(p.drumangle,1,size(r,1))')* pi/180));
d.pointing.xp_buddy = ...
    single(rb .* cos((thetab+repmat(p.drumangle,1,size(r,1))')* pi/180));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function maskground = get_groundmask(d,p,ind,expt,source_az,source_el)
  
% Ground + SPT + mast mask is all calculated using apparent az/el
% RELATIVE to beam centroid center.  Previously this used the 
% result of the 2D gaussian fits, but since we don't do that 
% step til later, let's trust the user-input source_az,source_el
% as good enough estimates to the beam center to make the mask.
% Definitely check this.
maskground = ones(size(d.pointing.az_ap),'uint8');
az = d.pointing.az_ap - source_az;
el = d.pointing.el_ap - source_el;

% Ground: mask out everything below el = -1.5
maskground(el < -1.5) = 0;

% SPT (Keck only): az between [-7,-2] and el < 1
if strcmp(expt,'keck')
  maskground(el < 1 & az > -7 & az < -2) = 0;
end

% Mast: az between [-0.3 0.3]
% How much to eat into the main beam in el depends on the beam size
% Semi-arbitrary, based on looking at the maps
% The mast is mostly obvious in Keck 220 maps - hard to see in B3
% 2016.  However let's be safe!  It has a similar beamwidth to Keck
% 220 so we'll go with that distance.
% Keck 270s -- match 220s for now, but may move mask closer in later.
for ii = 1:size(az,2)
  switch expt
    case 'keck'
      switch p.band(ii)
        case 100
          ellim = el(:,ii) < -1.2;
        case 150
          ellim = el(:,ii) < -1;
        case {210,220,270}
          ellim = el(:,ii) < -0.5;
      end
    case 'bicep3'
      ellim = el(:,ii) < -0.5;
  end
  azlim = az(:,ii) > -0.4 & az(:,ii) < 0.4;
  maskground(azlim & ellim,ii) = 0;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = rename_fields(d,p,expt,run)

% Find the chop reference
% Also keep utcfast time if B3/Keck
switch expt
  case 'bicep2'
    switch run
      case 'r81'
        b.ref = single(d.antenna0.pmac.fast_aux_input(:,2));
      otherwise
        b.ref = single(d.antenna0.pmac.fast_aux_input(:,4));
    end
  case 'keck'
    b.ref = single(d.antenna0.pmac.fast_aux_input(:,1)); % Keck
    b.t = d.antenna0.time.utcfast(:,1);
    b.t2 = d.antenna0.time.utcfast(:,2);
  case 'bicep3'
    b.t = d.antenna0.time.utcfast(:,1);
    b.t2 = d.antenna0.time.utcfast(:,2);
    switch run
      case 'hb3r5'
	b.ref = single(d.antenna0.pmac.fast_aux_input(:,1)); 
      case {'b3r6','b3r8'}
	b.ref = single(d.antenna0.pmac.fast_aux_input(:,4));
      case 'b3r9'
	if b.t(1) < date2mjd(2018) % 2017 campaign, only had chopper
          b.ref = single(d.antenna0.pmac.fast_aux_input(:,2));
	elseif b.t(1) < date2mjd(2019) % 2018 campaign w/ BSNS then chopper
	  b.ref = single(d.antenna0.pmac.fast_aux_input(:,1));
	end
    end
end  

% Detector feedback variable name
switch expt
  case 'keck'
    if ~isempty(d.mce0.data.fb)
      b.fb = double(d.mce0.data.fb);
    elseif ~isempty(d.mce0.data.raw)
      b.fb = double(d.mce0.data.raw);
    end
  case {'bicep2','bicep3'}
    if ~isempty(d.mce0.data.fb)
      b.fb = double(d.mce0.data.fb);
    elseif ~isempty(d.mce0.data.raw)
      b.fb = double(d.mce0.data.raw);
    end
end

% Keck 2013/14: sign flip for D2 in rx3
switch expt
  case 'keck'
    if b.t(1) > date2mjd(2013) & b.t(1) < date2mjd(2015)
      if rxNum == 3
	b.fb = -b.fb;
      else
	xx = find(p.rx == 3);
	b.fb(:,xx) = -b.fb(:,xx);
      end
    end
end

b.pointing = d.pointing;

return
