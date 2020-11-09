function tods = rps_load_tods(sch,schind,scans,p,rpsopt,dirname,chans,PW)
% function tods = rps_load_tods(sch,i,scans,p,rpsopt,dirname,chans,PW)
% Scans should be a 1xN array of scan numbers. Doesn't have to be
% sequential.

if size(scans,1) > 1
    scans = scans';
end

scans0 = scans;

if ~exist('PW','var')
    PW = false;
end

if PW
    if length(schind)==length(scans)
        PLOAD = true;
    else
        error('list of schedules and scans need to be the same length!')
    end
end

% if size(scans,2) > 1
%     tods = {};
% else
%     tods = [];
% end
tods = {};
for k = 1:length(schind)
    % Load in all timestreams for the input rasters in a schedule.
    i = schind(k);
    schname = sch{i}.name(end-15:end-4);
    
    
    if PW
        scans = scans0(k);
    end
    for j = scans
        fname = [dirname 'tods/tod_', schname, sprintf('_%03i.mat', j)];
        try
            load(fname);
            % Sometimes the tod is empty. Skip it.
            if isfield(rpstod,'az')
                if true %size(scans,2) > 1
                    tods{end+1} = rpstod; % Allows us to skip faulty timestreams if needed.
                else
                    tods = rpstod;
                end
            else
                fprintf('Scan %03i appears empty. Skipping...\n',j)
            end
        catch ME
            if strcmp(ME.identifier,'MATLAB:load:couldNotReadFile')
                fprintf('Could not load: %s\nSkipping...\n',fname)
            else
                disp(sprintf('Caught exception: %s\nWe''ll try to continue anyway.',ME.identifier))
            end
        end
        
        disp(['loaded ' fname '...'])
    end
    if length(tods)==0
        error('No data in TODs!\n')
    end
    
    
    % Convert to x/y focal plane coordinates by using the boresight parameters
    % r and theta which are identically zero.
    bs.r = 0;
    bs.theta = 0;
    bs.chi = 0;
    bs.chi_thetaref = 0;
    bs.drumangle = 0;
    %bs.drumangle = -unique(p.drumangle); %Change this so that this works for multi-drum exps...
    
    % Check if any specific channels are wanted.
    if exist('chans','var')
        channel = chans;
    else
        channel = tods{1}.ch;
        for i = 2:length(tods)
            channel = intersect(channel,tods{i}.ch);
        end
        
        % Find any channel among TOD's.
        for i = 2:length(tods)
            channel = union(channel,tods{i}.ch);
        end
        ch_count = zeros(size(channel));
        
        % Cut out channels that have less than 75% of rasters
        for i = 1:length(tods)
            ch_count = ch_count + ismember(channel,tods{i}.ch);
        end
        channel = channel(ch_count>=length(tods)*0.75);
        
    end
    
    % Add NaN timestreams to the channels that pass, but don't have
    % everything.
    %keyboard()
    for i = 1:length(tods)
        for j = 1:length(channel)
            if ~any(tods{i}.ch==channel(j))
                dummydata = NaN(size(tods{i}.az));
                tods{i}.ch(end+1) = channel(j);
                tods{i}.todcos(:,end+1) = NaN(size(tods{i}.todcos(:,1)));
                tods{i}.todsin(:,end+1) = NaN(size(tods{i}.todsin(:,1)));
                tods{i}.todquad(:,end+1) = NaN(size(tods{i}.todquad(:,1)));
            end
        end
    end
    
    % Check if TOD's have any channels at all.
    if isempty(channel)
        fprintf('Row %02i doesn''t have any channels. Skipping...\n',ii)
        return
    end
    %fprintf('Fitting a total of %i channels. This could take awhile.\n', length(channel))
    
    % rpstod.x and rpstod.y are actually det-centered coords (x' and y')
    % Change this to boresight-centered coordinates (x and y)
    
    rot = 1:length(tods);
    for j = 1:length(tods)
        [r, theta, psi] = ...
            keck_beam_map_pointing(tods{j}.az, tods{j}.el, tods{j}.dk, ...
            rpsopt.mount,rpsopt.mirror, rpsopt.source, bs);
        
        tods{j}.x = 2 * sind(r / 2) .* cosd(theta) * 180 / pi;
        tods{j}.y = 2 * sind(r / 2) .* sind(theta) * 180 / pi;
        tods{j}.phi = psi;
        %Grab source angles as well.
        rot(j)  = tods{j}.scan.abs_rot-rpsopt.grid.horizontal;
        tods{j}.rot = rot(j);
    end
    
end

