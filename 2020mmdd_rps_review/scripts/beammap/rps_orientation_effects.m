function parm_all = rps_orientation_effects(param,fd,p,C,PLOT)
% function rps_orientation_effects(a,b,offs,pol,C)


%% Investigate the systematic bias on psi due to RPS orientation

% Our model from Colin's posting
% A(theta) = [(cosd(psi)*cosd(b)+sind(psi)*sind(a)*sind(b))*cosd(th)+sind(psi)*cosd(a)*sind(th)].^2

simdirname  = 'rps_data/2019_mirror/ort_effects_sims/';
figdir = 'postings/2019mmdd_rps_dk_dep/figs3/';

% Start by grabbing the rough per-tile mean and std.
cutchans_a = unique(fd.ch(fd.inda==1));
cutchans_b = unique(fd.ch(fd.inda~=1));
if 0
ang = fd.phi_d;
dk = unique(fd.dk);
mn = NaN(20,2);
sd = mn;
for i = 1:2
    
    % Grab phi from each dk
    
    for k = 1:20
        
        [a, xpol, datarms] = deal(NaN(2640,4));
        ch = repmat((1:2640)',1,4);
        res = repmat(Inf,2640,4);
        for j = 1:length(dk)
            chind = find(fd.dk==dk(j) & p.tile(fd.ch)==k);
            ch(fd.ch(chind),j) = fd.ch(chind);
            a(fd.ch(chind),j) = ang(chind);
            res(fd.ch(chind),j) = fd.res_var(chind);
            datarms(fd.ch(chind),j) = fd.data_rms(chind);
            xpol(fd.ch(chind),j) = fd.aparam(chind,2);
        end
        
        
        if i==1
            ind = ismember(ch(:,1),cutchans_a);
        else
            ind = ismember(ch(:,1),cutchans_b);
            
        end
        
        s = nanstd(a(ind,:),0,2);
        m = nanmean(a(ind,:),2);
        [wm,wv,neff] = wmean(a(ind,:),1./res(ind,:),2);
        [wmrms] = wmean(datarms(ind,:),1./res(ind,:),2);
        ind2 = wm~=0 & wmrms ~=0;
        wm = wm(ind2);
        wmrms = wmrms(ind2);
        [wmwm, wvwm] = wmean(wm,1./wmrms);
        mn(k,i) = wmwm;
        sd(k,i) = sqrt(wvwm);
    end
end
end
% Grab estimated phi_d's by dk.
ang = fd.phi_d;
dk = unique(fd.dk);
for i = 1:2
    
    % Grab phi from each dk
    
    
    [a, xpol, datarms] = deal(NaN(2640,4));
    ch = repmat((1:2640)',1,4);
    res = repmat(Inf,2640,4);
    for j = 1:length(dk)
        chind = find(fd.dk==dk(j));
        ch(fd.ch(chind),j) = fd.ch(chind);
        a(fd.ch(chind),j) = ang(chind);
        res(fd.ch(chind),j) = fd.res_var(chind);
        datarms(fd.ch(chind),j) = fd.data_rms(chind);
        xpol(fd.ch(chind),j) = fd.aparam(chind,2);
    end
    
end

% Initialize parameters
th = -180:30:180;
sig = 0;
%sig2 = 0.1;
%offs = -2;
dk = [0 45 90 135]+1.25;
freepar.free = [1 1 1 1 1];
freepar.lb = [-360 -1 -1e4 -1e4 -1e4];
freepar.ub = [360 1 1e4 1e4 1e4];

parm_all = NaN(2640,4,5);
%parm_allx = [];

% Ignore the means from above for now. That might give us bias answers.
% Assume A/B pairs are exactly 90 degrees from each other.
mn = NaN(20,2);
sd = mn;

model = rps_get_ort_model(param,false);

for i = 1:20
    if unique(p.mce(p.tile==i)==0)
        mn(i,1) = normrnd(mce0a,sig);
        %mn(i,2) = atand(tand(mn(i,1)+89.5));
        mn(i,2) = normrnd(mce0b,sig);
    else
        mn(i,2) = normrnd(mce13b,sig);
        %mn(i,1) = atand(tand(mn(i,2)-89.5));
        mn(i,1) = normrnd(mce13a,sig);
    end
end

for j = 1:4
    dks = dk(j);
    %parm_psi = NaN(2640,1);
    if 0
        for i = [cutchans_a' cutchans_b']
            if strcmp(p.pol(i),'A')
                h=1;
            else
                h=2;
            end
            % Make sims
            
            k = p.tile(i);
            psi = dks+normrnd(mn(k,h),sd(k,h));
            A = ((cosd(psi)*cosd(el)+sind(psi)*sind(az)*sind(el))*cosd(th)+sind(psi)*cosd(az)*sind(th)).^2;
            guess = [psi 0.000 0 0 1/2];
            
            % Estimate parameters
            [aparam, aerr, agof, astat, acov] = matmin('chisq',...
                guess, freepar,'rps_get_mod_model',A,1,th);
            
            %parm_psi(i) = aparam(1)-dks;
            aparam(1) = aparam(1)-dks;
            
            parm_all(i,j,:) = aparam';
        end
    else
        for k = 1:20
            for h = 1:2
                psi = dks+mn(k,h);
                A = ((cosd(psi)*cosd(el)+sind(psi)*sind(az)*sind(el))*cosd(th)+sind(psi)*cosd(az)*sind(th)).^2;
                guess = [psi 0.000 0 0 1/2];
                
                % Estimate parameters
                [aparam, aerr, agof, astat, acov] = matmin('chisq',...
                    guess, freepar,'rps_get_mod_model',A,1,th);
                
                aparam(1) = aparam(1)-dks;
                if h==1
                ind = p.tile==k & strcmp(p.pol,'A');
                else
                ind = p.tile==k & strcmp(p.pol,'B');    
                end
                for i = find(ind)'
                parm_all(i,j,:) = aparam';
                end

            end
        end
    end
end

if 0
aznum = sign(az)*floor(abs(az));
azdec = floor(mod(az,1)*100);
elnum = sign(el)*floor(abs(el));
eldec = floor(mod(el,1)*100);

pname = [simdirname sprintf('orteff_sim_az_%i_%i_el_%i_%i_withoffset',aznum,azdec,elnum,eldec)];

simparam.az = az;
simparam.el = el;
simparam.fd = fd;
simparam.parm_all = parm_all;
save(pname,'simparam')

if 0
    
    aznum = sign(az)*floor(abs(az));
    azdec = floor(mod(az,1)*100);
    elnum = sign(el)*floor(abs(el));
    eldec = floor(mod(el,1)*100);
    
    pname = [simdirname sprintf('orteff_sim_az_%i_%i_el_%i_%i',aznum,azdec,elnum,eldec)];
    
    load(pname)
    parm_all = simparam.parm_all;
    
end

a1 = parm_all(:,1,1);
a2 = parm_all(:,2,1);
a3 = parm_all(:,3,1);
a4 = parm_all(:,4,1);
end


%% Show phi-fit as a function of DK
figdir = '~/postings/2019mmdd_rps_dk_dep/figs3/';
% best-fit so far:
%param = [4.6896 -4.111 88.803 -1.784 -0.905 88.314 -0.682 88.627 -0.982 88.458];
%pname = [figdir sprintf('orteff_sim_bestfit',aznum,azdec,elnum,eldec)];

% zero az/el offset
param = [2 -2.5 90-1.327 -1.784 -0.905 90-1.856 -0.682 88.627 -0.982 90-1.69];
pname = [figdir sprintf('orteff_sim_az_%i_%i_el_%i_%i',aznum,azdec,elnum,eldec)];


model = rps_get_ort_model(param,false);

az = param(1);
el = param(2);

ang.a0 = param(3);
ang.b0 = param(4);
ang.a1 = param(5);
ang.b1 = param(6);
ang.a2 = param(7);
ang.b2 = param(8);
ang.a3 = param(9);
ang.b3 = param(10);


fig = figure(1);
clf(fig)
set(fig,'Position',[1000 -100 [500 400]*2])
set(fig,'PaperPosition',[500 0 [5 4]*1.5])
dk = unique(fd.dk);

aznum = sign(az)*floor(abs(az));
azdec = floor(mod(az,1)*100);
elnum = sign(el)*floor(abs(el));
eldec = floor(mod(el,1)*100);


mk = {'--',''};
for j = 1:4
    subplot(4,1,j)
    hold off
    
    for i = 1:2
        if data(1,j,i)>10
            C = 90;
        else
            C = 0;
        end
        plot(dk',data(:,j,i)-C,['r' mk{i}])
        hold on
       
    end
    
    for i = 1:2
        if model(1,j,i)>10
            C = 90;
        else
            C = 0;
        end
        plot(dk',model(:,j,i)-C,mk{i})
        hold on
        
    end
    for i = 1:2
        if data(1,j,i)>10
            C = 90;
        else
            C = 0;
        end
        plot(dk',data(:,j,i)-C,['r' '.'])
        hold on
    end
    
    for i = 1:2
        if model(1,j,i)>10
            C = 90;
        else
            C = 0;
        end
        plot(dk',model(:,j,i)-C,'.')
        hold on
        
    end
    
    
    
    legend('PolA data','PolB data','PolA Sim','PolB Sim','Location','northeastoutside')
    title(sprintf('MCE%i Az=%i El=%i PolA=%2.2f PolB=%2.2f',j-1,aznum,elnum,...
        ang.(['a' sprintf('%i',j-1)]),ang.(['b' sprintf('%i',j-1)])))
    grid on
    xlim([0,140])
    ylim([-2 0])
    
    ylabel('phi_d (^o)')
end
xlabel('Dk (^o)')
%saveas(fig,pname,'png')

%% Plot only data.
% best-fit so far:


fig = figure(1);
clf(fig)
set(fig,'Position',[1000 -100 [500 400]*2])
set(fig,'PaperPosition',[500 0 [5 4]*1.5])
dk = unique(fd.dk);


pname = [figdir sprintf('orteff_real')];

mk = {'--',''};
for j = 1:4
    subplot(4,1,j)
    hold off
    
    for i = 1:2
        if data(1,j,i)>10
            C = 90;
        else
            C = 0;
        end
        plot(dk',data(:,j,i)-C,['r' mk{i}])
        hold on
        
    end
    
    for i = 1:2
        if data(1,j,i)>10
            C = 90;
        else
            C = 0;
        end
        plot(dk',data(:,j,i)-C,['r' '.'])
        hold on
        
    end
    
        legend('PolA','PolB','Location','northeastoutside')
        title(sprintf('MCE%i Weighted phi_d per DK',j-1))
    
    grid on
    xlim([0,140])
    ylim([-2 0])
    
    ylabel('phi_d (^o)')
end
xlabel('Dk (^o)')
saveas(fig,pname,'png')

%% Plot hist of best-fit residual


fig = figure(1);
clf(fig)
set(fig,'Position',[1000 0 [420 420]])
set(fig,'PaperPosition',[500 0 [4 4]])
pname = [figdir 'orteff_res_hist'];
hist(reshape(data,1,[])-reshape(model,1,[]))
xlim([-0.2 0.2])
grid on
xlabel('PolA/B weighted mean - best fit (^o)')
ylabel('N')
title('Best fit residual histogram')
saveas(fig,pname,'png')


%%

fig = figure(1);
clf(fig)
set(fig,'Position',[1000 0 [420 420]])
set(fig,'PaperPosition',[500 0 [4 4]])
pname = [figdir 'orteff_pairdiff_hist'];
z = [];
for j = 1:4
x = (data(:,j,1)+(data(:,j,2)-90))/2;
y = (ang.(['a' sprintf('%i',j-1)]) + (ang.(['b' sprintf('%i',j-1)])-90))/2;
z = [z; x-y];
end
hist(z)
xlim([0 0.3])
xlabel('Estimated bias in pair-diff angle per dk per mce')
ylabel('N')
grid on
saveas(fig,pname,'png')

