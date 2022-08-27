function posting_plots_rps_analysis_2()

clc; clear all;
%load('z:/dev/rps/rps_beam_fits_type11_rerun_cut.mat')
load('z:/dev/rps/rps_beam_fits_type5_rerun_cut.mat')
%load('z:/dev/rps/rps_beam_fits_rerun_all_cut.mat')
figdir = 'c:/Users/James/Documents/GitHub/postings/2022mmdd_rps_pointing/figs/';
load('z:/dev/rps/fpu_data_obs.mat')
load('z:/dev/rps/pm.mat')
load('z:/dev/rps/source_fit_data.mat')
inrange = @(A,B,C) B <= A & A <= C;
outrange = @(A,B,C) A <= B | C <= A;
addpath('z:/dev')
addpath('z:/pipeline/beammap')
addpath('z:/pipeline/util')
addpath('z:/dev/rps')


prx = 2*sind(p.r/2).*cosd(p.theta)*180/pi;
pry = 2*sind(p.r/2).*sind(p.theta)*180/pi;

clc
% Moon-derived mirror params -- Use only these and nothing else!
mirror = struct();
mirror.height = 1.4592;
mirror.tilt= 44.8870;
mirror.roll = -0.070;
rpsopt.mirror = mirror;

rpsopt.source.distance = 195.5;
% Fit for the source params given our mirror info:
source = rps_fit_source(fd,rpsopt,p,'');
rpsopt.source = source;


% With new mirror and source parameters, update the pointing.
[fd.x,fd.y,phi] = beam_map_pointing_model(fd.az_cen,fd.el_cen,fd.dk_cen,model,'bicep3',mirror,source,[]);
fd.x = reshape(fd.x,size(fd.ch));
fd.y = reshape(fd.y,size(fd.ch));
fd.resx = reshape(prx(fd.ch),size(fd.ch))-fd.x;
fd.resy = reshape(pry(fd.ch),size(fd.ch))-fd.y;
[resth, resr] = cart2pol(fd.resx,fd.resy);
[fd.resx_rot, fd.resy_rot] = pol2cart(resth-fd.dk_cen*pi/180,resr);

% Make a median subtraction array.
fd.phi_medsub = fd.phi;
for chind = 1:2640
    ind = fd.ch==chind;
    if ~isempty(find(ind))
        fd.phi_medsub(ind) = fd.phi_medsub(ind)-nanmedian(fd.phi(ind));
    end
end


%%
% Order the a/b phis based on schedule
[phis, xpols, phis_err, xs, ys] = deal(NaN(len,2640));
for schind = 1:len
    for chind = 1:2640
        ci = find(fd.ch==chind & ismember(fd.schnum,scheds{schind}));

        if ~isempty(ci)

            ch_chi = atand(tand(p.theta(chind)+p.chi(chind)));

            val = reshape(fd.phi(ci),[],1);
            xp = reshape(fd.xpol(ci),[],1);
            err = reshape(fd.phi_err(ci),[],1);

            %err(isnan(val))=1e10;
            xs(schind,chind) = mean(fd.x(ci));
            ys(schind,chind) = mean(fd.y(ci));
            phis(schind,chind) = mean(val);%wmean(val,1./err,1);
            xpols(schind,chind) = mean(xp);%wmean(xp,1./err,1);
            phis_err(schind,chind) = nanmin(err);

        end
    end
end



% Calculate pair-diff angles
% Loop over channels to account for MCE0;
[phi_pair, xpol_pair, phi_pair_err] = deal(NaN(len,2640));
for schind = 1:len
    for chind = 1:length(inda)

        pha = phis(schind,inda(chind));
        phb = phis(schind,indb(chind));
        ea = xpols(schind,inda(chind));
        eb = xpols(schind,indb(chind));

        phi_pair_err(schind,inda(chind)) = max(phis_err(schind,[inda(chind) indb(chind)]));
        if p.mce(inda(chind))~=0

            [phi_pair(schind,inda(chind)), xpol_pair(schind,inda(chind))] = calc_pair_diff_pol(pha,phb,ea,eb);
        else

            [phi_pair(schind,inda(chind)), xpol_pair(schind,inda(chind))] = calc_pair_diff_pol(phb,pha,eb,ea);
        end

        if 0%~inrange(phi_pair(schind,inda(chind))+2.5,-2.5,2.5)
            phi_pair(schind,inda(chind)) = NaN;
            xpol_pair(schind,inda(chind)) = NaN;
            phi_pair_err(schind,inda(chind)) = 1e10;
        end
    end
    %ind = abs(phi_pair(schind,:)-nanmean(phi_pair(schind,:)))<1;
end
