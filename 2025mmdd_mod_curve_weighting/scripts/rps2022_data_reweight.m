function fd = rps2022_data_reweight()

%%

load('dev/paper_data/rps_workspace_for_plotting.mat')
%%
%wtype = 'uniform';
wtype = 'oneoveramp';
%wtype = 'dadeoveramp';

%%
if 1
load(fullfile('dev','rps2022_data_for_BKXVIII_final'));

fd.parm_refit = {};
% Recycle parameter bounds
% [psi eps N1 N2 A]
freepar.free = [1 1 1 1 1];
freepar.lb = [-20 -1 -1e4 -1e4 0];
freepar.ub = [110 1 1e4 1e4 1e4];


tic
for chind = 1:length(fd.ch)


    % First guess.
    %guess = [fd.phi(chind) 0.003 0 0 max(fd.bparam{chind}(6:end))/2];
    switch wtype
        case 'uniform'
            new_weight = ones(size(fd.bparam{chind}(6:end)));
        case 'oneoveramp'
            new_weight = sqrt(0.005^2 + 0.03^2*cosd(fd.rot{chind}+fd.phi_s(chind)+fd.phi(chind)).^4);
        case 'dadeoveramp'
            dade = sqrt(sind(fd.rot{chind}+fd.phi_s(chind)+fd.phi(chind)).^2+1e-10);
            sig = sqrt(0.005^2 + 0.03^2*cosd(fd.rot{chind}+fd.phi_s(chind)+fd.phi(chind)).^4);
            new_weight = sig./dade;
    end
    
    % Estimate parameters
    [aparam, ~,~,~,~] = matmin('nanchisq',...
        [fd.phi(chind) 0.003 0 0 max(fd.bparam{chind}(6:end))/2],...
        freepar,...
        'rps_get_mod_model',...
        fd.bparam{chind}(6:end),...
        new_weight,...
        fd.rot{chind}+fd.phi_s(chind));

    fd.parm_refit{chind} = aparam;
    c = num2cell(aparam);
    [fd.phi(chind),fd.xpol(chind),fd.n1(chind),fd.n2(chind),fd.Amp(chind)] = deal(c{:});
end
toc
fd = get_pair_params(fd,ind0,ind90);
save(['dev/rps2022_data_refit_' wtype],'fd')

else
load(['dev/rps2022_data_refit_' wtype])

end

%%
fd_uniform = load('dev/rps2022_data_refit_uniform');
fd_uniform = fd_uniform.fd;
[fd_uniform, phis_uniform, phi_pair_uniform, xpols_uniform, poleffs_uniform,n1s_uniform,n2s_uniform,amps_uniform] = get_pol_params_per_obs(fd_uniform,p,scheds);


fd_refit = load(sprintf('dev/rps2022_data_refit_%s',wtype));
fd_refit = fd_refit.fd;
[fd_refit, phis_refit, phi_pair_refit, xpols_refit, poleffs_refit,n1s_refit,n2s_refit,amps_refit] = get_pol_params_per_obs(fd_refit,p,scheds);




%%
phi_pair_tilesub_refit = phi_pair_refit;
phi_pair_tilesub_uniform = phi_pair_uniform;

for tileind = 1:20
    idx = p.tile==tileind;% & ismember(1:2640,ind0)';
    phi_pair_tilesub_refit(:,idx) = phi_pair_refit(:,idx)-nanmedian(nanmedian(phi_pair_refit(:,idx),1),2);
    phi_pair_tilesub_uniform(:,idx) = phi_pair_uniform(:,idx)-nanmedian(nanmedian(phi_pair_uniform(:,idx),1),2);
    
end



% 

lims = {[-0.015, 0.025],[-3.25,-1],[-1,1]*0.25,[0.975,1.015]};
ylims = {[0,900],[0,120],[0,150],[0,250]};
vals = {'xpols','phi_pair','phi_pair_tilesub','poleffs'};
clabs = {'','[Deg]','[Deg]',''};
types = {'refit','uniform'};
fnames = {wtype,'uniform'};
samps = 25;
figdir = 'postings/2025mmdd_mod_curve_weighting/figs/';
printres = 300;
for validx = 3%1:length(vals)
    for tidx = 1:length(types)
        vname = [vals{validx} '_' types{tidx}];
        fname = [vals{validx} '_' fnames{tidx}];
        eval(['val = nanmedian(' vname ',1);'])

        eval(['val0 = ' vname ';']);
        
       
        
        fig = jfigure(5,[600,550]);
        clf()
        
        val(ind0) = nanmean([val(ind0);val(ind90)]);
        val(ind90) = val(ind0);
        fig.Position(1:2) = [10,369];
        plot_tiles(val,p,'clim',lims{validx},'title',strrep(vname,'_','\_'),'fig',5,'pair','mean');
        fname = sprintf('tile_plots_%s.png',fname);
        exportgraphics(fig,fullfile(figdir,fname),"Resolution",printres)
        
        vname = [vals{validx} '_' types{tidx}];
        fname = [vals{validx} '_' fnames{tidx}];
        eval(['val = nanmedian(' vname ',1);'])
        l = lims{validx};
        idx = val<=l(1) | val>=l(2);
        val(idx) = NaN;
        
        fig = jfigure(4,[400,400]);
        clf()
        fig.Position(1:2) = [10,369];
        edges = lims{validx}(1):diff(lims{validx})/samps:lims{validx}(2);
        N = histc(val,edges);
        bar(edges,N,'histc')
        xlim(lims{validx})
        ylim(ylims{validx})
        grid on
        
        title({strrep(vname,'_','\_'),sprintf('md: %0.3f | det2det std: %0.3f | : obs2obs std: %0.4f',nanmedian(val),nanstd(val),nanmedian(nanstd(val0,[],1)))})
        fname = sprintf('hists_%s.png',fname);
        exportgraphics(fig,fullfile(figdir,fname),"Resolution",printres)
    end
end

