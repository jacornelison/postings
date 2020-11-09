function bm_demod(bmopt)
% function bm_demod(bmopt)
% 
% Function to demodulate beam map timestreams.
% Data products are saved to a user-requested directory.
% The data product is bm_*.mat, which includes all of 
% of the demodulated timestreams.  There is one file saved for every
% corresponding arc file. 
%
% NOTE: AFTER BEAMMAP PIPELINE UPGRADE 07/2018, THIS CODE IS NOW
% OUTDATED.  USE BM_TIMESTREAM INSTEAD.
%       
% INPUTS (should all be passed in with bmopt)
%
%   t1:          start time (UTC), e.g. '2014-Mar-08:12:00:00'
%   t2:          stop time (UTC)
%   expt:        'bicep2','keck','bicep3'
%   run:         BICEP2 runs: 'r5','r8'
%                Keck runs:   'k0r4', 'k0r4t2' (tile 2 only), highbay
%                             'k083' = k7 run 2013-10, highbay
%                             'k201?' = Keck ffbm @ Pole (12,13,14,15)
%                             'k201?fsl' = Keck sidelobe @ Pole (13,14)
%                BICEP3 runs: 'hb3r5' = B3 run 5, highbay
%                             'b3r6' = B3 run 6 @ Pole, 2015-02
%                             'b3r8' = B3 run 8 @ Pole, 2016-02
%                             'b3r9' = B3 run 9 @ Pole, 2017-02
%   rxNum:       0 default, or specific Keck rx number, or 'all'
%   notlastscan: 0 default.  1 keeps the last scan if it's
%                incomplete, does not demodulate that scan, and
%                attaches it to the following file.  Must start from
%                the beginning of the beam map since the code is not
%                smart enough to figure out whether the beginning of
%                the file is a complete scan.  Slow! 
%   demodtype:   'square' (default), 'choplet'
%   demoddir:    directory in which to save demodulated data (normally
%                beammaps/demod)
%
% OPTIONAL
%
%   demodshift:  number of samples by which to shift ref signal - we have
%                hard-coded default values for various times, but can 
%                change this if we want.  Gets saved to demodded file.
%   pm:          pointing model parameters (overrides defaults)
%   v7p3:        forces saving in v7.3 (1 to use)
%   update:      option to not reprocess existing demodulated files. 
%                Default = 0 will re-process all files, 
%                = 1 only process missing files

% Parse bmopt and set defaults
t1 = bmopt.t1;
t2 = bmopt.t2;
expt = bmopt.expt;
run = bmopt.run;
if ~isfield(bmopt,'rxNum');
  rxNum = 0;
  else
  rxNum = bmopt.rxNum;
end
if ~isfield(bmopt,'notlastscan');
  notlastscan = 0;
else
  notlastscan = bmopt.notlastscan;
end
if ~isfield(bmopt,'demodtype');
  demodtype = 'square';
else
  demodtype = bmopt.demodtype;
end
if ~isfield(bmopt,'demoddir')
  demoddir = 'beammaps/demod';
else
  demoddir = bmopt.demoddir;
end
if ~isfield(bmopt,'demodshift')
  demodshift = [];
else
  demodshift = bmopt.demodshift;
end
if ~isfield(bmopt,'v7p3')
  bmopt.v7p3 = 0;
else
  bmopt.v7p3 = 1;
end
if ~isfield(bmopt,'update')
  update = 0;
else
  update = 1;
end


% Get a list of the relevant files - special case for highbay runs
switch run
  case {'k0r4','k0r4t2','k083','hb3r5'}
    files = list_arc_files('arc/',t1,t2,'.dat.gz');
  otherwise
    files = list_arc_files('arc/',t1,t2);
end 

p = get_array_info(files{1}(1:8));
  
% For concatenating previous files
e = [];

bmopt

% Loop over arc files:
for fl = 1:length(files)
  
  % if update, check if demodulated file already exists
  if update == 1
    flstr = [demoddir 'bm_' files{fl}(1,1:15) '_' demodtype '.mat']; 
    if exist_file(flstr) 
      disp(sprintf('arcfile %s already exists, skipping to next',flstr));
      continue
    end
  end
  
  tic
  load_file = char(strcat('arc/',files(fl)));
  disp(['Loading file: ' load_file]);
  % Load the data
  d = load_arc(load_file);
  disp(['File loaded (' num2str(toc) ' seconds)']);
  
  % Concatenate to previous file if requested
  if notlastscan
    dtmp = d;
    if ~isempty(e)
      clear d;    
      d = structcat(1,[e,dtmp]);
    end
    clear dtmp;
    dtmp = d;
  end

  % Data valid flag:
  mark = d.array.frame.features;

  % Check the data mode
  switch rxNum
    case 'all'
      dm = d.mce0.rc1.data_mode(logical(mark));
    otherwise
      dm = eval(sprintf('d.mce%d.rc1.data_mode(logical(mark))',rxNum));
  end
  
  % Check to see that there is valid data:
  if numel(dm) == 0
    disp(['File ' load_file ' is empty'])
    continue 
  end

  % Remove dropped samples
  syncant = d.antenna0.syncBox.sampleNumber;
  switch rxNum
    case 'all'
      syncmce = d.mce0.syncBox.sampleNumber(:,1);
    otherwise
      syncmce = eval(sprintf('d.mce%d.syncBox.sampleNumber(:,1)',rxNum));
  end
  
  % Identify a skip as any instance when they don't match
  d.skips = syncant ~= syncmce;

  % Generate field scans
  % Fast/slow ratio
  sampratio = ...
      length(d.antenna0.pmac.fast_az_pos)/length(d.antenna0.pmac.az_pos); 
  try
    % When is feature bit 1 on?
    fs = find_scans(d.antenna0.tracker.scan_off(:,1),...
	bitand(d.array.frame.features,2^1),sampratio,'keepend');
  catch
    warning('No field scans found in this file.  Moving to next file.')
    continue
  end
  disp(['Found scans (' num2str(toc) ' seconds)']);

  % Be sure there isn't a crazy turnaround right at the beginning where the
  % start of the scan happens after the beginning:
  fs = structcut(fs,fs.s < fs.e);
  % Also check that there are enough samples between the start/end of scan
  % Default Keck xfer function is 22 samples long, give a little buffer
  % Check for B3!
  fs = structcut(fs,fs.e - fs.s > 30);
  % If the above step killed the fs structure, go to next arcfile
  if isempty(fs.s)
    continue
  end
  
  % Keck/B3 options
  switch expt
    % Anything that isn't B2!
    case {'keck','bicep3'}
      for ii = 1:length(fs.e)
	if fs.e(ii) > length(d.mce0.cc.data_rate) 
	  fs.e(ii) = length(d.mce0.cc.data_rate); 
	end
      end
      % Keck pole-specific bad file
      if strcmp(files(fl),'20130224_214337.dat.gz')
	fs = structcut(fs,[1 2 3 5 6 7]);
      end
      % Keck at Pole/B3 has multiple MCEs - concatenate these
      switch run
	case {'k2012','k2013','k2013fsl','k2014','k2014fsl','k2015','k2016','k2017','k2018'} 
	  for m = {'mce1','mce2','mce3','mce4'}
	    if(isfield(d,m))
	      d.mce0 = structcat(2,[d.mce0,getfield(d,m{1})]);
	      d = rmfield(d,m);
	    end
	  end
	case {'hb3r5','b3r6','b3r8','b3r9'} % Same for B3
	  for m = {'mce1','mce2','mce3'}
	    if(isfield(d,m))
	      d.mce0 = structcat(2,[d.mce0,getfield(d,m{1})]);
	      d = rmfield(d,m);
	    end
	  end
      end % Concatenate switch
      % Deconvolve if at Pole - something weird with B3 for now
      switch run
	case {'k2012','k2013','k2013fsl','k2014','k2014fsl','k2015','k2016','k2017','k2018','b3r8','b3r9'}
	  d = get_xferfunc(d,fs,0,0); 
	  d = deconv_scans(d,p,1);
	  disp(['Deconv scans done (' num2str(toc) ' seconds)']);
      end
  end % Keck/B3 options
  
  d.t = d.antenna0.time.utcfast(:,1);
  
  % Generate pointing model
  % This isn't with the input parsing section since it wants the "files" var
  if ~isfield(bmopt,'pm');
    mjd = date2mjd(str2num(files{1}(1:4)),str2num(files{1}(5:6)),...
	str2num(files{1}(7:8)));
    bmopt.pm =  get_pointing_model(mjd,0,d,fs);
  end
  bmopt.pm.az_tilt_ha = 0;
  bmopt.pm.az_tilt_lat = 0;
  bmopt.pm.el_tilt = 0;
  
  pm = bmopt.pm;

  % Inv pointing to get hor coords.  
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
  disp(['Pointing model done (' num2str(toc) ' seconds)']);
  
  % Get the data rate
  data_rate = mode(double(d.antenna0.syncBox.data_rate));

  % Far sidelobe stuff
  switch run
    case {'k2013fsl','k2014fsl'}
      d = cut_the_crap(d,fs,rxNum);
  end 

  % Keep only fields pos, fb, skips, and ref chop 
  % 2015-03-06, change pos to pointing
  d = rename_fields(d,rxNum,run,p,expt);

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
  
  % CORE DEMODULATION ROUTINE
  disp('Demodulation start');
  d = demod_blocks(d,fs,sqw,dm,run,demodtype,expt,demodshift);
  disp(['Demodulation done (' num2str(toc) ' seconds)']);
  
  % Cut out dropped samples if we're not B3 at the high bay!
  switch run
    case {'hb3r5'}
      disp('Not removing skips')
    otherwise
     % d = structcut(d,~d.skips);
  end
  
  % Zero out transients:
  %trans = diff([d.ind; 0]);
  %d.sin(trans > 100,:) = 0;
  %d.cos(trans > 100,:) = 0;
  
  % Save out the data:
  % If there is a user-input demod shift, save it in the filename
  if ~isempty(demodshift)
    demodstr = ['_shift' int2str(demodshift)];
  else
    demodstr = [];
  end
  switch expt
    case 'bicep2'
      svstr = [demoddir '/bm_' files{fl}(1,1:15) '_' demodtype demodstr]; 
      save(svstr,'d')
    case 'keck'
      switch rxNum
	case 'all'
	  svstr = [demoddir '/bm_' files{fl}(1,1:15) '_' demodtype demodstr];
	  if bmopt.v7p3
	    save(svstr,'d','-v7.3');
	  else
	    save(svstr,'d');
	  end
	otherwise
	  svstr = [demoddir '/bm_' files{fl}(1,1:15) '_rx' ...
		int2str(rxNum) '_' demodtype demodstr];
	  save(svstr,'d');
      end % rxNum switch
    case 'bicep3'
      svstr = [demoddir '/bm_' files{fl}(1,1:15) '_' demodtype demodstr];
      if bmopt.v7p3
	save(svstr,'d','-v7.3');
      else
	save(svstr,'d')
      end
  end
  system(['chmod g+w ' svstr '.mat']);
  disp(['File ' svstr ' saved (' num2str(toc) ' seconds)']);
end

return % Main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = demod_blocks(d,fs,sqw,dm,run,demodtype,expt,demodshift)

% Get rid of NaNs in timestream and sqw for Keck
% How necessary is this???
switch expt
  case 'keck'
    indices = find(isnan(d.fb(:,300))); % Is this just a random channel?
    for ii = length(indices):-1:1
      d.fb(indices(ii),:) = [];
      d.ref(indices(ii)) = [];
      %d.pos(indices(ii),:) = []; % Phased out for d.pointing
      d.pointing.hor.az(indices(ii)) = [];
      d.pointing.hor.el(indices(ii)) = [];
      d.pointing.hor.dk(indices(ii)) = [];
      %d.com(indices(ii),:) = [];
      d.skips(indices(ii)) =[];
      sqw(indices(ii)) = [];
    end
end

switch demodtype
  case 'choplet'
    d.mce0.data.fb = d.fb;
    % What circshift do we want?
    % 2012-02
    if d.t(1) < date2mjd(2013) 
      shiftme = 2;
    % 2013-02
    elseif d.t(1) < date2mjd(2013,2,14) 
      shiftme = 7;
    % 2013-02 to 2013-03
    elseif d.t(1) > date2mjd(2013,2,18) & d.t(1) < date2mjd(2013,3,10) 
      shiftme = 2;
    end
    % Apply the shift
    if ~isempty(demodshift)
      shiftme = demodshift;
    end
    d.antenna0.pmac.fast_aux_input(:,1) = circshift(d.ref,[shiftme,0]);    
    %d.antenna0.time.utcfast=cat(2,d.t,d.t2);
    d.t = utc2datenum(cat(2,d.t,d.t2));
    d = demod_choplet(d);
    %d.cos = d.mce0.data.fb;
    d = rmfield(d,{'antenna0','fb','mce0','t','t2'})  
    
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
    
    % Shift sqw then clean chop ref before demodding
    d.sqw = shiftchopref(sqw,shiftme);
    ch = cleanchopref(d.sqw);
    [e.cos,e.sin,e.ind] = square_demod(d.fb,ch.refclean_hires,'highres');
    
    if length(e.sin(:,1)) > length(e.cos(:,1))
      e.sin = e.sin(2:length(e.cos(:,1)) + 1,:);
    elseif length(e.sin(:,1)) < length(e.cos(:,1))
      e.sin = [e.sin; ...
	zeros(length(e.cos(:,1)) - length(e.sin(:,1)),size(e.sin,2))];
    end
    d = rmfield(d,{'t','t2'});
    
    % Downselect indices that are within the field scans
    % CHECK VS KNOWN GOOD KECK
    cc = [];
    for ii = 1:length(fs.sf)
      keep = e.ind(e.ind > fs.sf(ii) & e.ind < fs.ef(ii));
      cc = [cc; keep];
    end
    
    % This doesn't do the downsampling correctly since we added the
    % 'pointing' struct to d.  Ugly fix for now...
    pointing = d.pointing;
    d = rmfield(d,'pointing');
    % Downsample
    d = structcut(d,cc);
    pointing.hor = structcut(pointing.hor,cc);
    % And put pointing back
    d.pointing = pointing;
    
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
    
    d.cos = e.cos; 
    d.sin = e.sin; 
    d.ind = e.ind;
    d.demodshift = shiftme;
    d = rmfield(d,'sqw');

end % Standard demod

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = rename_fields(d,rxNum,run,p,expt)

% Find the chop reference
switch expt
  case 'bicep2'
    switch run
      case 'r81'
	b.ref = double(d.antenna0.pmac.fast_aux_input(:,2)); % BICEP2
      otherwise
	b.ref = double(d.antenna0.pmac.fast_aux_input(:,4)); % BICEP2 
    end
  case 'keck'
    b.ref = double(d.antenna0.pmac.fast_aux_input(:,1)); % Keck
    % Keep time
    b.t = d.t;
    b.t2 = d.antenna0.time.utcfast(:,2);
  case 'bicep3'
    switch run
      case 'hb3r5'
	b.ref = double(d.antenna0.pmac.fast_aux_input(:,1)); 
      case {'b3r6','b3r8'}
	b.ref = double(d.antenna0.pmac.fast_aux_input(:,4));
      case 'b3r9'
	if d.t(1) < date2mjd(2018) % 2017 campaign, only had chopper
          b.ref = double(d.antenna0.pmac.fast_aux_input(:,2));
	elseif d.t(1) < date2mjd(2019) % 2018 campaign w/ BSNS then chopper
	  b.ref = double(d.antenna0.pmac.fast_aux_input(:,1));
	end
    end
    b.t = d.t;
    b.t2 = d.antenna0.time.utcfast(:,2);
end

% Pointing info

%b.pos(:,1) = d.pointing.hor.az;
%b.pos(:,2) = d.pointing.hor.el;
%b.pos(:,3) = d.pointing.hor.dk;

%b.com(:,1) = d.command.az;
%b.com(:,2) = d.command.el;
%b.com(:,3) = d.command.dk;

b.pointing = d.pointing;
b.skips = d.skips;

% Assign variable names
switch expt
  case 'keck'
    switch rxNum
      case 'all'
	if ~isempty(d.mce0.data.fb)
	  b.fb = double(d.mce0.data.fb);
	elseif ~isempty(d.mce0.data.raw)
	  b.fb = double(d.mce0.data.raw);
	end
      otherwise
	xx = find(p.rx == rxNum);
	if ~isempty(d.mce0.data.fb(:,xx))
	  b.fb = double(d.mce0.data.fb(:,xx));
	elseif ~isempty(d.mce0.data.raw(:,xx))
	  b.fb = double(d.mce0.data.raw(:,xx));
	end
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
    if d.t(1) > date2mjd(2013) & d.t(1) < date2mjd(2015)
      if rxNum == 3
	b.fb = -b.fb;
      else
	xx = find(p.rx == 3);
	b.fb(:,xx) = -b.fb(:,xx);
      end
    end
end

% Keep source position if we're observing the moon:
if isfield(d,'source_az')
  b.source_az = d.source_az;
  b.source_el = d.source_el;
  b.source_pa = d.source_pa;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sqrw = deglitch_sqw(sqw)
glitch = sqw;
deb = diff([0; sqw]);
for jj = 2:length(sqw)-1
  if deb(jj) == -1 && deb(jj+1) == 1
    glitch(jj) = 1;
  end
  if deb(jj) == 1 && deb(jj+1) == -1
    glitch(jj) = 0;
  end
end
sqrw = glitch;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pm = get_command_coords(t,expt)
% Can someone explain why this is called pmi???
pmi = get_pointing_model(t);

% John changed the encoder_zeros in schedlib on Jan 17 2011 based on star
% pointing results.  These values are hard-coded into schedlib - I can't
% figure out a way to get them from the data except to hard code it here
% ...you mean hard-coded into pointing.init, right?  Cause that's where it is

switch expt
  case 'bicep2'
    if t > date2mjd(2011,01,17)
      pm.az_zero = -322.493;
      pm.el_zero = -82.436;
    else
      pm.az_zero = -322.79;
      pm.el_zero = -82.398;
    end
    pm.el_tilt = 0;
    pm.az_tilt_lat = 0;
    pm.az_tilt_ha = 0;
  case 'keck'
    pm.az_zero = 1.175575;
    pm.el_zero = 0.106166;
    pm.el_tilt = -0.000030;
    pm.az_tilt_lat = -0.006825;
    pm.az_tilt_ha = 0;
  case 'bicep3' % In pointing.init, 2015-01-07 SAH
    pm.az_zero = -6.8;
    pm.el_zero = -82.398;
    pm.el_tilt = 0;
    pm.az_tilt_lat = 0;
    pm.az_tilt_ha = 0;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = cut_the_crap(d,fs,rxNum)

for i = 1:length(fs.s)
  s = fs.sf(i); 
  e = fs.ef(i);
  fbstd = nanstd(d.mce0.data.fb(s:e,:));
  d.mce0.data.fb(s:e,fbstd > 2000) = NaN;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
