function mparam = rps_fit_mirror_params(sch,p,rpsopt,ch)

% How do we want to do this?
% Find a good a/b pair from each tile, trying to get an even distribution
% on the focal plane.
% Load pertinent timestreams
% Fit for beam centers in az/el space
% Fit mirror params
%   Feed az_0/el_0/dk into kbmp.m
%   Use observed beam centers for p
%   Use residual x_0/y_0 for chi-square minimization

% Loop something like
% bs.r = 0;
%     bs.theta = 0;
%     bs.chi = 0;
%     bs.chi_thetaref = 0;
%     bs.drumangle = 0;
%     mount.aperture_offr = 0;

% [r th psi] = keck_beam_map_pointing(az_0,el_0,dk,...
%       rpsopt.mount,rpsopt.mirror,rpsopt.source,bs);

% Grab az/el beam centers
load('man_cuts')
good_chans = find(man_cuts);

prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;
dirname = 'rps_data/2018/';
fd = [];
fd.ch = [];
fd.dk = [];
fd.azel = [];

% Manually cut tiles so far:
% [1 2 3 4 8 13 11]
%rot = -180:30:180;
for i = 1:length(sch)
    schname = sch{i}.name(end-15:end-4);
    for j = 1:sch{i}.nrows
        rpsparam_file = [dirname, 'params/param_', schname, sprintf('_%02i.mat', j)];
        try
            load(rpsparam_file)
            if ~isempty(i_param) & ~isfield(i_param{1},'map_param')
                for k = 1:length(i_param)
                    %if any(p_ind.rgl100==i_param{k}.ch)
                    if any(good_chans==i_param{k}.ch) & any( [1:20] == p.tile(i_param{k}.ch))
                        fd.azel = [fd.azel; i_param{k}.azel_cen];
                        
                        fd.ch = [fd.ch; i_param{k}.ch];
                        fd.dk = [fd.dk; i_param{k}.dk];
                        
                    end
                end
            end
        catch
            disp(['Something went wrong. Skipping: ', rpsparam_file])
            continue
        end
    end
end


% cut channels
if ~exist('ch','var') | isempty(ch)
    ch = fd.ch;
    azel = fd.azel;
    dk = fd.dk
else
    ind = ismember(fd.ch,ch);
    ch = fd.ch(ind);
    azel = fd.azel(ind,:);
    dk = fd.dk(ind);
end

% first guess
guess = [rpsopt.mirror.tilt rpsopt.mirror.roll];
% Set up boundaries

freepar.free = ones(size(guess));
freepar.lb = [30 -20];
freepar.ub = [80 20];
% Estimate parameters
% Our data is the observed beam centers, prx and pry
[param, err, gof, stat, cov] = matmin('chisq',...
    guess, freepar,	'rps_get_mirror_model',[prx(ch); pry(ch)],[],azel,dk,rpsopt);

mparam.param = param;
mparam.err = err;
mparam.gof = gof;
mparam.stat = stat;
mparam.cov = cov;
mparam.xy = reshape(rps_get_mirror_model(param,azel,dk,rpsopt),[],2);
mparam.ch = ch;
mparam.dk = dk;



