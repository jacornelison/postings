function  out = color_constraint_b1xb2(varargin)
% Uses B1(100) and B1(100)xB2(150) to limit a foreground contribution
% to B2(150) BB.  Works with 0751x1651 simset.
%
% 09Jan2014 JPF Use H-L likelihood, treating 150 and 100 as if T, E
% 11Feb2014 JPF New B1xB2 sims
% 19Feb2014 JPF Major reorganization, many options
% 05Mar2014 JPF Try MCMC marginalization
% 09Mar2014 JPF Back to profiling by default (faster, same answer),
%               now more self-contained.  Paper plot defaults.
% 15Mar2014 JPF Allow same constraint for EE spectrum
% 16May2014 JPF Remove spectra plotting.  Add sub_lensing option.  Fix CMB scaling in xDust.
% 20May2014 JPF Add color correction for realistic spectra.

%% ARGUMENT AND ENVIRONMENT PREP

p = inputParser;
p.addParamValue('bins',2:6); % bin(s) to plot
p.addParamValue('n',-4:0.05:4); % assumed spectral indices
p.addParamValue('plot',0); % should I make plots?
p.addParamValue('n_offdiagonal',1); % max ell bin separation in cov matrix (makes little difference)
p.addParamValue('marginalize',0); % marginalize over nuisance parameters (instead of profile; if MCMC is off, this is a product of 1D marginalizations)
p.addParamValue('MCMC',0); % Marginalize with a multi-dimensional MCMC
p.addParamValue('prior',[-4,4]); % Range for Bayesian integral
p.addParamValue('sim_r',1); % fiducial model with nonzero r?
p.addParamValue('hack_dof',0); % Hack B1 to have same degrees of freedom as B2? CAREFUL!!!
p.addParamValue('hack_depth',0); % Hack B1 to have more depth? CAREFUL!!!
p.addParamValue('sym_bayes',0); % Symmetric (rather than shortest) credible interval
p.addParamValue('show_95',0); % Show 95% confidence lines in plots?
p.addParamValue('like_ratio_title',1); % Show likelihood ratio values in title
p.addParamValue('EE',0); % Use EE instead of BB
p.addParamValue('sub_lcdm',1); % Subtract lensing from BB, LCDM from EE
p.addParamValue('use_spectra',0); % Use measured band pass spectra
p.addParamValue('is_odyssey',1); % Use default paths
p.parse(varargin{:});

% move a few parameters into the main workspace for convenience
n = p.Results.n; % spectral indices to plot
nbins = length(p.Results.bins);
if nbins==1
    bintxt = ['Bin ' int2str(p.Results.bins(1))];
else
    bintxt = ['Bins ' int2str(p.Results.bins(1)) ' - ' int2str(p.Results.bins(end))];
end

%% SPECTRAL INDEX PREP
% Fundamental constants
k = 1.38e-23; % J/K (Boltzmann's constant)
h = 6.63e-34; % J-s (Planck's constant)
c = 299793e3; % m/s (speed of light)
T_cmb = 2.7255; % Kelvin

nuB1 = 96.0e9; % Hz
nuB2 = 149.8e9; % Hz
BW = 0.25; % receiver fractional bandwidth

% Band pass spectra in convenient structure form
specB1 = [];
specB2 = [];
if p.Results.use_spectra
    % use real bands (for now captured from B1/B2 paper plots)
    datB1 = importdata('/Users/Jeff/matlab/BICEP2/0751x1651/BICEP1_100.txt');
    datB2 = importdata('/Users/Jeff/matlab/BICEP2/0751x1651/BICEP2_150.txt');
    datB2.data(:,1) = 60+(datB2.data(:,1)-50)*(250-60)/(250-50);% hack for bad digitization!
    specB1.spec = @(nu) interp1(datB1.data(:,1)*1e9,datB1.data(:,2),nu,[],0);
    specB2.spec = @(nu) interp1(datB2.data(:,1)*1e9,datB2.data(:,2),nu,[],0);
    specB1.domain = [min(datB1.data(:,1)) max(datB1.data(:,1))]*1e9;
    specB2.domain = [min(datB2.data(:,1)) max(datB2.data(:,1))]*1e9;
else
    % assume tophat bands
    specB1.domain = nuB1*[(1-BW/2) (1+BW/2)];
    specB2.domain = nuB2*[(1-BW/2) (1+BW/2)];
    specB1.spec = @(nu) inrange(nu,specB1.domain(1),specB1.domain(2)); 
    specB2.spec = @(nu) inrange(nu,specB2.domain(1),specB2.domain(2)); 
end

% BLACKBODY RELATIONS (one polarization, 100% efficiency)
% Spectral radiance (power per unit area, solid angle, and frequency)
Iplanck = @(v,T) (h/c^2)*v.^3./(exp(h*v./(k*T))-1);
%Assume a single-moded receiver (throughput A-Omega = lambda^2)
% Power absorbed by a single-moded receiver with band sp
Pplanck = @(sp,T) quad(@(v) sp.spec(v).*Iplanck(v,T).*(c./v).^2, sp.domain(1),sp.domain(2));
% Differential power near T_cmb for conversions to Kcmb (ignore normalization)
dPcmb = @(sp) (Pplanck(sp,T_cmb+0.001) - Pplanck(sp,T_cmb-0.001))/0.002;


% FOREGROUND POWER LAW RELATIONS (ignore constant of proportionality
% See arXiv:1207.3675 (PSM) and 1312.5446 (Planck dust) for notation.
% My n is their beta: the power law of the (single-moded) T_antenna.
% So T(v) ~ v^n  ~ lambda^2 I(v)  [single-moded receiver!]
% and I(v) ~ v^(n+2) [spectral radiance]
% SYNCH: n ~ -3
% DUST: greybody (beta~1.5-1.65, T~20K) looks like I~v^(3.4), so n~1.4
Ifg = @(v,n) (v/150e9).^(n+2); % [normalized for convenience]
% Assume a single-moded receiver (throughput A-Omega = lambda^2)
% That adds an additional factor of v^(-2) to the integrand over the
% band, and that integral adds a further factor of v for fractional BW.
% (which I do analytically to avoid matrix dimension annoyances)
%Pfg = @(v0,n) integral(@(v) Ifg(v,n).*(c./v).^2, v0*(1-BW/2),v0*(1+BW/2));
%Pfg = @(v0,n) integral(@(v) (v/150e9).*n, v0*(1-BW/2),v0*(1+BW/2));
%Pfg = @(v0,n) (v0/150e9).^(n+1); % ignore normalization constant
Pfg = @(sp,n) quad(@(v) sp.spec(v).*Ifg(v,n).*(c./v).^2, sp.domain(1),sp.domain(2));

% equivalent index for "standard" Planck dust: T=19.6K, beta=1.6 (Aumont talk)
Igrey = @(v,beta,Td) Iplanck(v,Td).*(v/150e9).^beta;
% Assume a single-moded receiver (throughput A-Omega = lambda^2)
% That adds an additional factor of v^(-2) to the integrand
Pgrey = @(sp,beta,Td) quad(@(v) sp.spec(v).*Igrey(v,beta,Td).*(c./v).^2,...
    sp.domain(1),sp.domain(2));


%% LOAD DATA
% convenience indices: spectra
TT=1; TE=2; EE=3; BB=4; TB=5; EB=6;
% convenience indices: experiments and colors
B2=1; B1=2; B1150=3; X=4; X150=5; XB1=6;
% reduc_final output
% from Clem's final/0751x1651
if p.Results.is_odyssey
    datdir = '/n/bicepfs1/bicep2/pipeline/final/0751x1651/';
else
    datdir = '/Users/Jeff/matlab/BICEP2/0751x1651/';
    %fnameFiducial = '/Users/Jeff/matlab/BICEP2/0751x1651/camb_planck2013_r0_lensing.fits';
end
if p.Results.hack_depth
    %fname = 'real_a_filtp3_weight3_gs_dp1102_jack0_real_a_filtp3_weight3_gs_dp1000_jack0_pureB_matrix_0p3.mat';
    %fname = 'real_a_filtp3_weight3_gs_dp1102_jack0_real_a_filtp3_weight3_gs_dp1000_jack0_pureB_matrix_0p1.mat';
    fname = 'real_a_filtp3_weight3_gs_dp1102_jack0_real_a_filtp3_weight3_gs_dp1000_jack0_pureB_matrix.mat';
else
    fname = 'real_a_filtp3_weight3_gs_dp1102_jack0_real_a_filtp3_weight3_gs_dp1000_jack0_pureB_matrix.mat';
end
d = load(fullfile(datdir,fname));
if p.Results.hack_depth
    nsim = size(d.r(B1).noi,3);
    rlz = randi(nsim); % pick a random realization to work with
    d.r(B1).real(:,BB) = d.r(B2).real(:,BB) + squeeze(d.r(B1).noi(:,BB,rlz)) - d.r(B1).db(:,BB);
    d.r(X).real(:,BB) = d.r(B2).real(:,BB) + squeeze(d.r(X).noi(:,BB,rlz)) - d.r(X).db(:,BB);
end
% factors of amplification for various power laws (x = P96/P150)
x_n = @(n) (Pfg(specB1,n)./dPcmb(specB1)) ./ (Pfg(specB2,n)./dPcmb(specB2));
n_cmb = fzero(@(n) x_n(n)-1,0); % effective power law of the CMB across these frequencies
disp(['CMB equivalent index: n = ' num2str(n_cmb,'%2.2f')]);

%betaD=1.6; TD=19.6;
%betaD=1.62; TD=19.7; % Planck XXXI (2013) all-sky intensity
betaD=1.52; TD=18.7; % Planck XXII (2014) galactic intensity
betaDP=1.63; TD=18.7; % Planck XXII (2014) galactic polarization
% xDust = Pgrey(nuB1,betaD,TD)/Pgrey(nuB2,betaD,TD);  % OLD: power ratio
xDust = (Pgrey(specB1,betaD,TD)/dPcmb(specB1))/(Pgrey(specB2,betaD,TD)/dPcmb(specB2)); % NEW: CMB ratio
n_dust = fzero(@(n) x_n(n)-xDust,1.5); % effective power law of fiducial dust across these frequencies
disp(['Dust (beta=' num2str(betaD,'%2.2f') ', TD=' num2str(TD,'%3.1f') ' K) equivalent index: n = ' num2str(n_dust,'%2.2f')]);
xDustP = (Pgrey(specB1,betaDP,TD)/dPcmb(specB1))/(Pgrey(specB2,betaDP,TD)/dPcmb(specB2)); % NEW: CMB ratio
n_dustP = fzero(@(n) x_n(n)-xDustP,1.5); % effective power law of fiducial dust across these frequencies
disp(['DustPol (beta=' num2str(betaDP,'%2.2f') ', TD=' num2str(TD,'%3.1f') ' K) equivalent index: n = ' num2str(n_dustP,'%2.2f')]);

xset = zeros(size(n)); % explicit loop, since not a vectorized function
for ii = 1:length(xset)
    xset(ii) = x_n(n(ii)); % multipliers for power spectra
end

if p.Results.EE
    spidx = EE;
    splabel = 'EE';
else
    spidx = BB;
    splabel = 'BB';
end
if p.Results.sub_lcdm
    % hacked computation of lensing or SM spectrum to subtract it
    % Note sim has dimensions (bin, spectrum, realization)
    offset = squeeze(mean(d.r(B2).sim(p.Results.bins,spidx,:),3));
    % Load lensed-LCDM from official FITS file
    %fiducialModel = fitsread(fnameFiducial,'asciitable');
else
    offset = 0;
end

% Extract appropriate spectra: real, noise sims, s+n sims
ell = d.r(B2).l;
BB_B1 = squeeze(d.r(B1).real(:,spidx));
BB_B2 = squeeze(d.r(B2).real(:,spidx));
BB_X = squeeze(d.r(X).real(:,spidx));
snsetB1 = squeeze(d.r(B1).sim(:,spidx,:));
snsetB2 = squeeze(d.r(B2).sim(:,spidx,:));
snsetX = squeeze(d.r(X).sim(:,spidx,:));

%% ====== NOTES ======
% I want to compute likelihood as a function of n: L(n)
% From there I can use a chi2 test (-2 log L) or an integral
% For profile likelihoods this is easy: a fast optimization for each n.
% The integral scheme demands a multi-dim integral for each n (painful).
% In practice each n demands a ~5D integral (5-bin power spectrum)
% We can approach this with a MCMC, but even that is slow: a big chain for
% every n!
%
% Instead, just generate one long chain over all parameters and then histogram
% over n.  That uses our computational time better: we want to spend most
% of our time generating points near the optimum in n, not equal numbers at
% all n.  BUT this is a bit tricky: the likelihood function you get
% will depend on the binning.  Binning depends implicitly on having a
% measure, and thus on having a prior!  Binning in x vs. n changes
% your tail shape, for example, and associated likelihood ratios.
%
% So in practice it's probably best to throw one MC for each value of
% x.  BUT this ends up coming out essentially identical to the profile
% likelihood, which is much faster to compute!  So use that.

%% ================ MULTI BIN CALCULATIONS ================

if p.Results.marginalize
    if p.Results.MCMC
        disp('Multi-bin marginalization with MCMC');
        % prepare for H-L likelihood calculation
        res = hl_create_data_products_2BB(d.r([B2,B1,X]), p.Results.bins, p.Results.n_offdiagonal,p.Results.sim_r,p.Results.EE);
        % Hack degrees of freedom, if desired (DON'T USE THIS USUALLY!!)
        if p.Results.hack_dof % FIX THIS!!!!
            dofB2 = dof(1,:);
            dofB1 = dof(2,:);
            % hack B1 noise to approximately match B2 d.o.f.
            for kk=1:nbins
                res.M_f(3*(kk-1)+2,3*(kk-1)+2) = res.M_f(3*(kk-1)+2,3*(kk-1)+2) * dofB1(kk)/dofB2(kk); % B1
                res.M_f(3*(kk-1)+3,3*(kk-1)+3) = res.M_f(3*(kk-1)+3,3*(kk-1)+3) * sqrt(dofB1(kk)/dofB2(kk)); % Cross
            end
        end
        prep = hamimeche_lewis_prepare(res.C_fl, res.M_f);
        % NOTE: everything here is fixed EXCEPT res.C_l_hat
	if 0
	  % one MCMC for everyone
	  fun = @(C) hl_logL_atX(prep,res,C(2:end),C(1),1);
	  % Throw N random band power triplets following this distribution
	  N=1e5; Nchain = 8;
	  samp = zeros(N,Nchain,nbins+1);
	  %wid = [0.4 0.04*ones(1,nbins)];
	  wid = [2 0.02*ones(1,nbins)]; % choose width to be just larger than typical distribution width
	  parfor cc = 1:Nchain
	    tic;
	    x0 = [1; d.r(B2).real(p.Results.bins,spidx)] + ...
		randn(nbins+1,1).*wid'/4; % initial guess
	    samp(:,cc,:) = ...
		slicesample(x0,N,'logpdf',fun,'burnin',max(200,floor(N/10)),'thin',3,'width',wid);
	    toc;
	  end
	  samp = reshape(samp,[],nbins+1);
	  % form into a flattened estimate of L(n)
	  mlike = hist(samp(:,1),xset(end:-1:1)); % need binning to be increasing
	  mlike = log(mlike(end:-1:1)); % now corresponds to n
	  save hlMC.mat n xset mlike samp p
	else
	  % one MCMC per x value
	  mlike = zeros(length(xset),1); % preallocate
	  wid = 0.02;
	  N=3e4;
	  parfor cc = 1:length(xset)
	    disp(['Chain ' int2str(cc)  ' of ' ...
		  int2str(length(xset))]);
	    tic;
	    fun = @(C) hl_logL_atX(prep,res,C,xset(cc),1);
	    x0 = d.r(B2).real(p.Results.bins,spidx) + ...
		randn(nbins,1)*wid/10; % initial guess
	    % generate the MCMC
	    samp = ...
		slicesample(x0,N,'logpdf',fun,'burnin',max(200,floor(N/10)),'thin',1,'width',wid);
	    % compute the summed likelihood
	    temp = 0;
	    for dd=1:N
	      temp = temp + exp(fun(samp(dd,:)));
	    end
	    mlike(cc) = log(temp);
	    toc
	  end
	end
    else
        disp('Multi-bin marginalization not applied');
        disp('Multiplying single-bin marginalizations');
        % LOOP OVER BINS, FILL BIG MATRIX WITH LOG(L(x))
        allbins = p.Results.bins;
        mlike = zeros([length(xset),length(allbins)]); % preallocate
        for idx = 1:length(allbins)
            thisbin = allbins(idx);
            disp(['Bin ' int2str(thisbin)]);
            
            % prepare for H-L likelihood calculation
            res = hl_create_data_products_2BB(d.r([B2,B1,X]), thisbin, p.Results.n_offdiagonal,p.Results.sim_r,p.Results.EE);
            % Hack degrees of freedom, if desired (DON'T USE THIS USUALLY!!)
            if p.Results.hack_dof
                dofB2 = 2*res.C_fl(1,1)^2/res.M_f(1,1);
                dofB1 = 2*res.C_fl(2,2)^2/res.M_f(2,2);
                % hack B1 noise to approximately match B2 d.o.f.
                res.M_f(2,2) = res.M_f(2,2) * dofB1/dofB2;
                res.M_f(3,3) = res.M_f(3,3) * sqrt(dofB1/dofB2);
            end
            prep = hamimeche_lewis_prepare(res.C_fl, res.M_f);
            
            % Evaluate at desired spectral indices
            for ii=1:length(xset)
                fun = @(A) arrayfun(@(C) hl_logL_atX(prep,res,C,xset(ii),0),A);
                mlike(ii,idx) = log(quad(fun,0,1));
            end
            
        end
        
        % Combine all bins (product of likelihoods = sum of log-likelihoods)
        mlike = sum(mlike,2);
    end

else
    % PROFILE LIKELIHOOD
    % prepare for H-L likelihood calculation
    res = hl_create_data_products_2BB(d.r([B2,B1,X]), p.Results.bins, p.Results.n_offdiagonal,p.Results.sim_r,p.Results.EE);

    % compute H-L equivalent degree-of-freedom statistics
    % means
    cvec = reshape(res.C_fl,4,[]);
    cvec = cvec([1,4,2],:); % (B2, B1, cross) x (bins)
    % variances
    mvec = reshape(diag(res.M_f),3,[]); % (B2, B1, cross) x (bins)
    % dof (Immanuel's prescription)
    dof = 2*(cvec.^2)./mvec;
   
    % Hack degrees of freedom, if desired (DON'T USE THIS USUALLY!!)
    if p.Results.hack_dof % FIX THIS!!!!
        dofB2 = dof(1,:);
        dofB1 = dof(2,:);
        % hack B1 noise to approximately match B2 d.o.f.
        for kk=1:nbins
            res.M_f(3*(kk-1)+2,3*(kk-1)+2) = res.M_f(3*(kk-1)+2,3*(kk-1)+2) * dofB1(kk)/dofB2(kk); % B1
            res.M_f(3*(kk-1)+3,3*(kk-1)+3) = res.M_f(3*(kk-1)+3,3*(kk-1)+3) * sqrt(dofB1(kk)/dofB2(kk)); % Cross
        end
    end
    prep = hamimeche_lewis_prepare(res.C_fl, res.M_f);
    
    
    % preallocate
    mlike = zeros([length(xset),1]); % preallocate
    spec0 = d.r(B2).real(p.Results.bins,spidx); % starting guess is B2 spectrum
    for ii=1:length(xset)
        % Solve for best-fit spectrum at this x (e.g. profile likelihood)
        [specmax, logLval] = fminsearch(@(spec) -hl_logL_atX(prep,res,spec,xset(ii),1,offset),spec0);
        mlike(ii) = -logLval;
    end
end

% Profile likelihood
plike = -2*(mlike - max(mlike));
[MLE,int68,int95] = get_frequentist_interval(n,mlike);
sigH = int68(2)-MLE;
sigL = MLE-int68(1);
disp(['Profile estimate: ' num2str(MLE,'%3.2f') ...
    ' (+' num2str(sigH,'%3.2f') '/-' num2str(sigL,'%3.2f') ')']);
disp(['   95% confidence: [' num2str(int95,'%3.2f') ']']);

% PLOTS
if p.Results.plot
    % Profile likelihood
    figure;
    plot(n,plike,'k-','LineWidth',3);
    hold on;
    set(gca,'FontSize',16);
    grid();
    if p.Results.EE
        axis([-4,4,0,100]);
    else
        axis([-4,4,0,7]);
    end
    h1=plot(n_cmb*[1,1],ylim,'-','Color',[0,0.7,0],'LineWidth',4); % CMB
    h2=plot(int68(1)*[1,1],ylim,'b--','LineWidth',3); % 68%
    plot(int68(2)*[1,1],ylim,'b--','LineWidth',3);
    if p.Results.show_95
        h3=plot(int95(1)*[1,1],ylim,'m:','LineWidth',3); % 95%
        plot(int95(2)*[1,1],ylim,'m:','LineWidth',3);
    end
    xlabel('Spectral index (\beta)','FontSize',18);
    ylabel('-2 ln \lambda_{profile}','FontSize',16);
    title(['H-L, ' bintxt ': ' num2str(int95(1),'%3.2f') ' < \beta < ' ...
        num2str(int95(2),'%3.2f')],'FontSize',18);
    if p.Results.show_95
        legend([h1,h2,h3],{'CMB','68% CL','95% CL'},'Location','SouthEast');
    else
        legend([h1,h2],{'CMB','68% CL'},'Location','SouthEast');
    end
    text(2.1,3.1,[num2str(MLE,'%3.2f') ...
        '^{+' num2str(sigH,'%3.2f') '}_{-' ...
        num2str(sigL,'%3.2f') '}'],'FontSize',24);
end

% Compute Bayesian limits with flat prior
% Allowed interval is smallest containing 95% of total integral
% Select range of integration (flat prior)
prior = n>p.Results.prior(1) & n<p.Results.prior(2);
L = exp(mlike);
n0 = n(prior);
L0 = L(prior);
[MPE,int68B,int95B]=get_bayesian_interval(n,squeeze(mlike),prior(:), ...
    p.Results.sym_bayes);
sigH = int68B(2)-MPE;
sigL = MPE-int68B(1);
disp(['Bayesian estimate: ' num2str(MPE,'%3.2f') ...
    ' (+' num2str(sigH,'%3.2f') '/-' num2str(sigL,'%3.2f') ')']);
disp(['   95% credible: [' num2str(int95B,'%3.2f') ']']);

if p.Results.plot
    figure;
    plot(n0,L0/max(L0),'k-','LineWidth',3);
    hold on;
    set(gca,'FontSize',16);
    grid();
    axis([-4,4,0,1]);
    h1=plot(n_cmb*[1,1],ylim,'-','Color',[0,0.7,0],'LineWidth',4); % CMB
    h2=plot(int68B(1)*[1,1],ylim,'b--','LineWidth',3); % 68%
    plot(int68B(2)*[1,1],ylim,'b--','LineWidth',3);
    if p.Results.show_95
        h3=plot(int95B(1)*[1,1],ylim,'m:','LineWidth',3); % 95%
        plot(int95B(2)*[1,1],ylim,'m:','LineWidth',3);
    end
    xlabel('Spectral index (\beta)','FontSize',18);
    ylabel('Likelihood','FontSize',16);
    if p.Results.like_ratio_title
        title(['H-L, ' bintxt ': ' num2str(int95(1),'%3.2f') ' < \beta < ' ...
            num2str(int95(2),'%3.2f')],'FontSize',18);
    else
        title(['H-L, ' bintxt ': ' num2str(int95B(1),'%3.2f') ' < \beta < ' ...
            num2str(int95B(2),'%3.2f')],'FontSize',18);
    end
    if p.Results.show_95
        legend([h1,h2,h3],{'CMB','68% CL','95% CL'},'Location','NorthEast');
    else
        legend([h1,h2],{'CMB','68% CL'},'Location','NorthEast');
    end
    text(1.5,0.4,[num2str(MPE,'%3.2f') ...
        '^{+' num2str(sigH,'%3.2f') '}_{-' ...
        num2str(sigL,'%3.2f') '}'],'FontSize',24);
end

%% ASSEMBLE OUTPUTS FOR PLOTTING
out.opt = p.Results; % options for plotting
out.n = n;
out.n_cmb = n_cmb;
out.mlike = mlike;
out.plike = plike;
out.MLE = MLE;
out.int68 = int68;
out.int95 = int95;
out.MPE = MPE;
out.int68B = int68B;
out.int95B = int95B;
out.xset = xset;

end


% ========================================================

%% --- Helper functions ---

function logL = hl_logL_atX(prep,res,C_l,x,islog,offset)
% Evaluates log(L) or Lin H_L approx for spectrum C_l (n_ell vector) and
% scale factor x.
if nargin<6
    offset=0;
end
if nargin<5
    islog=1;
end
n_ell = length(C_l);
temp = zeros(2,2,n_ell);
temp(1,1,:) = C_l+offset;
temp(2,2,:) = x^2*C_l+offset;
temp(1,2,:) = x*C_l+offset;
temp(2,1,:) = temp(1,2,:);
logL = hamimeche_lewis_likelihood(prep,reshape(res.C_l_hat,2,2,[]), ...
    temp+res.N_l,islog);
end


function [MLE,int68,int95] = get_frequentist_interval(n,logL)
% Return the peak, 1-sigma interval, and 95% interval from a log-likelihood
% in 1 dimension using a -2*logL ~ chi2 approximation

stat = -2*logL;
% find MLE
[minstat,idxMLE] = min(stat);
MLE = n(idxMLE);
stat = stat - minstat; % Distributed as a chi2 with 1 d.o.f.

int68 = get_crossings(n,stat,1.0,MLE); % 1-sigma ~ 68%
int95 = get_crossings(n,stat,chi2inv(0.95,1),MLE); % 95%

end


function int = get_crossings(x,y,thresh,x0)
% get threshold crossings to the left and right of x0 (if they exist)
fun = @(xx) interp1(x,y-thresh,xx);
try
    int(1) = fzero(fun,[x(1),x0]);
catch except
    int(1) = NaN;
end
try
    int(2) = fzero(fun,[x0,x(end)]);
catch except
    int(2) = NaN;
end

end

function [MPE,int68,int95] = get_bayesian_interval(n,logL,prior,sym_interval)
% Return the peak, 1-sigma interval, and 95% interval from a log-likelihood
% in 1 dimension using a bayesian integral and a flat prior

if nargin<4
    sym_interval = 0;
end
if nargin<3
    prior = ones(size(logL));
end
L = exp(logL).*prior; % posterior
% find MPE
[maxL,idxMPE] = max(L);
MPE = n(idxMPE);

if sym_interval
    % find intervals with equal integrals in both tails
    % This need not contain the MPE for short interval widths
    Lcdf = cumsum(L)/sum(L); % cumulative distribution function
    [Lcdf0,idx,dummy] = unique(Lcdf,'last'); % make single-valued
    n0 = n(idx);
    fun = @(p) interp1(Lcdf0,n0,p); % inverse CDF
    val68 = chi2cdf(1.0,1);
    int68 = fun([0.5-val68/2, 0.5+val68/2]);
    int95 = fun([0.5-0.95/2, 0.5+0.95/2]);

else
    % find cutoff L's such that points of larger L contribute x% of integral
    % This is the highest-posterior interval, and should be the shortest
    % This interval must contain the MPE for any interval width
    [Lsort,Lidx] = sort(L);
    Llim68 = max(Lsort(find(cumsum(Lsort) < (1-chi2cdf(1.0,1))*sum(Lsort)))); % Likelihood of limit
    int68 = get_crossings(n,L,Llim68,MPE);
    Llim95 = max(Lsort(find(cumsum(Lsort) < (1-0.95)*sum(Lsort)))); % Likelihood of limit
    int95 = get_crossings(n,L,Llim95,MPE);
end

end
function result = hl_create_data_products_2BB(r, ellbins, n_offdiagonal,simr,EE)
%
%Inputs
%  r - 3-element results structure from reduc_final_comb
%  ellbins - ell bins to include in the likelihood
%  n_offdiagonal - if non-negative, only use this many offdiagonal elements of bandpower covariance matrix
%  simr - if true, use sim with r nonzero (r.simr, rather than r.sim)
%  EE - if true, use EE instead of BB
%
%Like hl_create_data_products, BUT BB-only from multiple maps
%
%Outputs
%  result.N_l - [n_field, n_field, n_ell] noise bias
%  result.C_fl - fiducial model bandpowers + noise bias
%  result.C_l_hat - real data bandpowers + noise bias
%  result.M_f - bandpower covariance matrix

  if nargin<5
      EE=0;
  end
  if nargin<4
      simr=0;
  end

  if EE
      fields='E';
      spidx=3;
  else
      fields='B';
      spidx=4;
  end

  n_ell = length(ellbins);
  n_sim = size(r(1).sim, 3);
  n_fields = length(fields);
  n_spectra = n_fields*(n_fields+1)/2; %power spectra types

  % Now account for multiple maps; just like multiple fields
  % Check number of dimensions in r is a triangular number
  % Adds a leading index to all data structures
  n_maps = (sqrt(8*length(r)+1)-1)/2;
  if floor(n_maps)~=n_maps
      % Not triangular!
      crash_not_implemented;
  end
  n_fields2 = n_fields*n_maps;
  n_spectra2 = n_fields2*(n_fields2+1)/2;
  
  %Ordering in r: 150, 100, cross
  %Spectra ordering in pipeline: 'TT','TE','EE','BB','TB','EB','ET','BT','BE'
  
  %Preallocate matrices
  N_l = zeros(n_fields2,n_fields2,n_ell);
  C_fl = zeros(n_fields2,n_fields2,n_ell);
  C_l_hat = zeros(n_fields2,n_fields2,n_ell);
  
  %Fill
  N_l(1,1,:) = r(1).db(ellbins, spidx); % 150x150
  N_l(2,2,:) = r(2).db(ellbins, spidx); % 100x100
  N_l(1,2,:) = r(3).db(ellbins, spidx); % 150x100
  N_l(2,1,:) = N_l(1,2,:);
  
  if simr
      C_fl(1,1,:) = mean(r(1).simr(ellbins, spidx, :), 3);
      C_fl(2,2,:) = mean(r(2).simr(ellbins, spidx, :), 3);
      C_fl(1,2,:) = mean(r(3).simr(ellbins, spidx, :), 3);
      C_fl(2,1,:) = C_fl(1,2,:);
  else
      C_fl(1,1,:) = mean(r(1).sim(ellbins, spidx, :), 3);
      C_fl(2,2,:) = mean(r(2).sim(ellbins, spidx, :), 3);
      C_fl(1,2,:) = mean(r(3).sim(ellbins, spidx, :), 3);
      C_fl(2,1,:) = C_fl(1,2,:);
  end
  
  C_l_hat(1,1,:) = r(1).real(ellbins, spidx);
  C_l_hat(2,2,:) = r(2).real(ellbins, spidx);
  C_l_hat(1,2,:) = r(3).real(ellbins, spidx);
  C_l_hat(2,1,:) = C_l_hat(1,2,:);
  
  %Signal+noise sims have had noise bias (and E->B) subtracted so add it back
  C_fl = C_fl + N_l;

  %Real data had noise bias subtracted so add it back
  C_l_hat = C_l_hat + N_l;

  %Covariance is unchanged by adding a constant so it doesn't matter whether we add noise bias before calculating the covariance matrix
  %Interleave the spectra in the order required for covariance matrix
  %Covariance matrix order is 150, 100, cross at each ell
  sims = zeros(n_ell * 3, n_sim);
  if simr
      sims(1:3:(n_ell .* 3), :) = r(1).simr(ellbins, spidx, :);
      sims(2:3:(n_ell .* 3), :) = r(2).simr(ellbins, spidx, :);
      sims(3:3:(n_ell .* 3), :) = r(3).simr(ellbins, spidx, :);
  else
      sims(1:3:(n_ell .* 3), :) = r(1).sim(ellbins, spidx, :);
      sims(2:3:(n_ell .* 3), :) = r(2).sim(ellbins, spidx, :);
      sims(3:3:(n_ell .* 3), :) = r(3).sim(ellbins, spidx, :);
  end
  
  M_f = cov(sims');
 
  %Trim offdiagonal elements of covariance matrix  
  if (n_offdiagonal >=0)
    M_f_trimmed = zeros(size(M_f));
    for iSpectrum = 1:n_spectra2
      for jSpectrum = 1:n_spectra2
        %Get each submatrix for a pair of spectrum types
        submatrix = M_f(iSpectrum:n_spectra2:(n_ell.*n_spectra2), ...
          jSpectrum:n_spectra2:(n_ell.*n_spectra2));
        %This is an n_ell * n_ell matrix so trim in ell
        submatrix = triu(tril(submatrix, n_offdiagonal), -n_offdiagonal);
        %Insert matrix back
        M_f_trimmed(iSpectrum:n_spectra2:(n_ell.*n_spectra2), ...
          jSpectrum:n_spectra2:(n_ell.*n_spectra2)) = submatrix;
      end
    end
    M_f = M_f_trimmed;
%    M_f = triu(tril(M_f, n_offdiagonal), -n_offdiagonal);
  end

  %Add systematic uncertainty after trimming covariance matrix.
  %The trimming is because we cannot estimate offdiagonal elements well with limited MC
  %However, the systematic uncertainty covariance is not limited by MC
  %so don't trim systematic uncertainty -- IDB 2013-04-08
  %Add systematic uncertainty
  if isfield(r, 'abscal_uncer') | isfield(r, 'beam_uncer')
    model = make_hl_vector(r.real, fields, ellbins, n_ell);
  end

  if isfield(r, 'abscal_uncer')
    gain_uncer = make_hl_vector(r.abscal_uncer, fields, ellbins, n_ell);
    M_f = M_f + (cvec(gain_uncer)*rvec(gain_uncer)).* (cvec(model)*rvec(model));
  end
  if isfield(r, 'beam_uncer')
    beam_uncer = make_hl_vector(r.beam_uncer, fields, ellbins, n_ell);
    M_f = M_f + (cvec(beam_uncer)*rvec(beam_uncer)).* (cvec(model)*rvec(model));
  end

  result.N_l = N_l;
  result.C_fl = C_fl;
  result.C_l_hat = C_l_hat;
  result.M_f = M_f;
end

function vec = make_hl_vector(source, fields, ellbins, n_ell)

    switch fields
      case 'T'
        vec = source(ellbins, 1);
      case 'E'
        vec = source(ellbins, 3);
      case 'B'
        vec = source(ellbins, 4);
      case 'EB'
        vec = zeros(n_ell .* 3, 1); 
      %Interleave the spectra in the order required for covariance matrix
      %Covariance matrix order is EE, BB, EB at each ell
        vec(1:3:(n_ell .* 3)) = source(ellbins, 3); %EE 
        vec(2:3:(n_ell .* 3)) = source(ellbins, 4); %BB
        vec(3:3:(n_ell .* 3)) = source(ellbins, 6); %EB
      case 'TB'
        vec = zeros(n_ell .* 3, 1); 
      %Interleave the spectra in the order required for covariance matrix
      %Covariance matrix order is TT, BB, TB at each ell
        vec(1:3:(n_ell .* 3)) = source(ellbins, 1); %TT 
        vec(2:3:(n_ell .* 3)) = source(ellbins, 4); %BB
        vec(3:3:(n_ell .* 3)) = source(ellbins, 5); %TB
      case 'TE'
        vec = zeros(n_ell .* 3, 1); 
      %Interleave the spectra in the order required for covariance matrix
      %Covariance matrix order is TT, EE, TE at each ell
        vec(1:3:(n_ell .* 3)) = source(ellbins, 1); %TT 
        vec(2:3:(n_ell .* 3)) = source(ellbins, 3); %EE
        vec(3:3:(n_ell .* 3)) = source(ellbins, 2); %TE

      case 'TEB'
      %Interleave the spectra in the order required for covariance matrix
      %Covariance matrix order is TT, EE, BB, TE, EB, TB at each ell
        vec = zeros(n_ell .* 6, 1); 
        vec(1:6:(n_ell .* 6)) = source(ellbins, 1); %TT 
        vec(2:6:(n_ell .* 6)) = source(ellbins, 3); %EE
        vec(3:6:(n_ell .* 6)) = source(ellbins, 4); %BB
        vec(4:6:(n_ell .* 6)) = source(ellbins, 2); %TE
        vec(5:6:(n_ell .* 6)) = source(ellbins, 6); %EB
        vec(6:6:(n_ell .* 6)) = source(ellbins, 5); %TB


      otherwise
        crash_not_implemented
    end
end
