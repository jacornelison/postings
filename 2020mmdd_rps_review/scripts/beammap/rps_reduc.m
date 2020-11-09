function rps_reduc(sch,i,j,p,p_ind,rpsopt,dirname)
% rps_reduc(sch,i,j,p,rpsopt,dirname)
%
% Wrapper function that runs rps_read and rps_fit_beam
% [Input]
%   
%
%
%
%
%
%
%
%
% rps_reduc(sch,i,j,p,rpsopt,dirname)

%Initialize paths and variables
schname = sch{i}.name(end-15:end-4);
coverage_file = [dirname, 'coverage_plots/coverage_', schname, sprintf('_%03i.eps', j)];
%disp(coverage_file)
rpstod_file = [dirname, 'tods/tod_', schname, sprintf('_%03i.mat', j)];
%disp(rpstod_file)

if isempty(p)
    mjdstr = mjd2datestr(sch{i}.scans(j).t1);
    datevect = [mjdstr(10:11), mjdstr(6:8), mjdstr(1:4)]; %Have to swap year and day because datestr is stupid.
    
    fprintf('.\n.\n.\nGetting array info\n.\n.\n.\n')
    arrinfodate = datestr(datevect , 'yyyymmdd');
    p = get_array_info(arrinfodate);
    
end

%% Find on-source channels.
if strcmp(p.expt,'keck') | strcmp(p.expt,'bicep3')
    fprintf('.\n.\n.\nFinding on-source channels\n.\n.\n.\n')
    [chans, min_dist] = rps_find_channels(sch{i}.scans(j), p, rpsopt, coverage_file);
elseif isfield(rpsopt, 'ch')
    chans = rpsopt.ch;
else
    % Do all channels
    chans = 1:length(p.gcp);
end

chans = intersect(chans, p_ind.rgl);
chans = sort(chans);

fprintf('Found %4.0f channels.\n',size(chans,2))
fprintf(['Saving coverage plot to:\n' coverage_file '\n'])


%% RPS REDUC
% Get TOD's
fprintf('.\n.\n.\nRunning RPS_read\n.\n.\n.\n')
rpstod = rps_read(sch{i}.scans(j), chans, p, rpsopt);
rpstod.min_dist = min_dist;
fprintf('.\n.\n.\nrps_read complete!\n.\n.\n.\n')
fprintf(['Saving reduc data to:\n' rpstod_file '\n.\n.\n.\n'])

%%%% SAVE
save(rpstod_file, 'rpstod', '-v7.3'); % -v7.3 compresses .mat file.

disp('Reduc complete!')
end
