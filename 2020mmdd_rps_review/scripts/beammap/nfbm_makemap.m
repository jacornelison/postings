function [px x y sq] = nfbm_makemap(savename, t1, t2, steps, type)
% function [px x y sq] = nfbm_makemap(savename, t1, t2, steps, type)
%
% Demodulates and filters raw time streams and returns maps
% 
% Inputs:
%    savename = name of .mat file where maps are saved
%    t1 = start of scan
%    t2 = end of scan

%    steps = number of steps (assumed to be the same for x and y)
%    type = option to select sin, cos, or quadrature sum returned from
%           demodulation.
%
% Old inputs:
%    x_enc_min = encoder count for beginning of slew in x
%    x_enc_max = encoder count for end of slew in x
%    y_enc_min = encoder count for beginning of slew in y
%    y_enc_max = encoder count for end of slew in y
%
% Outputs:
%    px(row,col,tile) = struct with elements 
%                       px.a_map = A polmap, px.b_map = B pol map
%    x = x coordinates in encoder counts
%    y = y coordinates in encoder counts
%    sq(gcp_ind) = structure with element sq.map = dark squid maps 

%-----------------------------------------
% EXAMPLES by JET
% K3_2:
% >>[px x y sq]=near_field_bm('2010-dec-29:09:53:30','2010-dec-29:14:14:16', 
%                             34849, 1038, -186148, 221034, 56, 'cos') 
% K1_4
% >>[px x y sq]=near_field_bm('2011-jan-22:15:05:15','2011-jan-22:19:27:50', 
%                             -112796, 107200,-113349, 110649, 56, 'quad')
%
% time from: >>t=(utc2time(b.antenna0.time.utcfast(2.7e4, :))/60/60-9)*60 
% y_encmin from: >>min(b.antenna0.pmac.fast_enc_pos(:,7))
% x_encmin from: >>min(b.antenna0.pmac.fast_enc_pos(:,8))
%
% After this run:
% >>nf_beam_plotter_K3_2(x,y,px,sq, 10 ,14, 'gauss', 'rotate', 3)
% >>nf_beam_plotter_K3_2(x,y,px,sq, 20 ,14, 'gauss', 'rotate', 3)
% >>nf_beam_plotter_v2(x,y,px,sq, 25 ,14, 'gauss')
%----------------------------------------------

% Load up all arcfiles between start and stop times
% Specified file extension so list_arc_files doesn't choke on .part files
files = list_arc_files('arc/',t1,t2,'.dat.gz');

% Filter parameters:
[bf af] = butter(1, 1/50, 'high');

% Loop through all arcfiles
for fl = 1:length(files)
   
  load_file = char(strcat('arc/', files(fl)))

  % We want feedback, sample number, clock, x/y encoder positions, chop ref
  % [6,7] are counted started from 0!!!!    
  b = load_arc(load_file,[],[],{'mce0.data.fb','antenna0.syncBox.sampleNumber','mce0.header.clock_counter','antenna0.pmac.fast_enc_pos[6,7]','antenna0.pmac.fast_aux_input[0]'});

  % Assign var names:
  fb = b.mce0.data.fb;
    
  % Shift indices:  GCP index 1:527 = matlab index 1:527, 
  %                 GCP index 0 = matlab index 528.
  fb = circshift(fb,[0,-1]);
    
  % Turn chop reference signal to a logical square wave
  % Threshold should be roughly half of the TTL that PMAC sees
  sqw = b.antenna0.pmac.fast_aux_input < 8000;
    
  % Translate encoder counts to x/y positions
  % This uses the conversion factor encoder counts -> inches (Velmex?)
  y = double(b.antenna0.pmac.fast_enc_pos(:,1)) / 16000;
  x = double(b.antenna0.pmac.fast_enc_pos(:,2)) / 16000;

  % Reform data to eliminate dropped bins
  ind_nz = b.mce0.header.clock_counter ~= 0 & ...
      b.antenna0.syncBox.sampleNumber ~= 0 | (y ~= 0 & x ~= 0);
  fb = fb(ind_nz,:);
  x = x(ind_nz);
  y = y(ind_nz);
  sqw = sqw(ind_nz);

  % Filter the fb:
  hp_a = filtfilt(bf, af, fb);
    
  % Reform the flags to be a clean chop, then figure out where the
  % transitions are:
  [cos_data,sin_data,edges] = square_demod(hp_a,sqw);

  if length(sin_data) > length(cos_data) 
    sin_data = sin_data(1:end-1,:);  
  end
  if length(sin_data) > length(cos_data) 
    sin_data = sin_data(1:end-1,:);  
  end

  quad = sqrt(cos_data.^2+sin_data.^2);
    
  part(fl).x_pos = x(edges)';
  part(fl).y_pos = y(edges)';
    
  if strcmp(type,'sin')
    part(fl).amp = sin_data';
    disp('using sin')
  elseif strcmp(type,'cos')
    part(fl).amp = cos_data';
    disp('using cos')
  else
    part(fl).amp = quad';
    disp('using quadrature sum')
  end
    
  plot(quad(:,12))
end

% Concatenate over all files
x_pos = [part(:).x_pos];
y_pos = [part(:).y_pos];
amp = [part(:).amp];
amp = amp';

% Make x/y grids
x = linspace(min(x_pos),max(x_pos),steps);
y = linspace(min(y_pos),max(y_pos),steps);

% Loop over tiles, rows, columns
for ti = 1:4
  for ro = 1:8
    for co = 1:8
      a_ind = det2gcp(ro,co,'A',ti);
      b_ind = det2gcp(ro,co,'B',ti);
      if a_ind == 0
	a_ind = 528;
      end
      % Grid up the map for each pixel
      px(ro,co,ti).a_map = grid_map(x_pos,y_pos,amp(:,a_ind),x,y);
      px(ro,co,ti).b_map = grid_map(x_pos,y_pos,amp(:,b_ind),x,y);
    end
  end
end

% List of dark squids -- only for 150 GHz, change for 100!
dsq = [1,34,67,100,133,166,199,232,265,298,331,364,397,430,463,496];
sq = [];
for ii = 1:length(dsq)
  sq(ii).map = grid_map(x_pos,y_pos,amp(:,dsq(ii)),x,y);
end

save(savename,'px','x','y','sq')

return


