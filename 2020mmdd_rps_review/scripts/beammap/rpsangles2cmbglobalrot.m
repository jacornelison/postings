function rpsangles2cmbglobalrot(phiQ,fitsfile,cutparamfile,finalfile)
% function rpsangles2cmbglobalrot(phiQ,finalfile,cutparamfile)
% For a given input of channels with measures polarization angles,
% calculate a back-of-the-envelope estimate of the global rotation angle
% without running map-based sims.


if isstr(fitsfile)
    disp('Loading final file...')
    almsT=read_fits_alms(fitsfile,1);
    almsE=read_fits_alms(fitsfile,2);
    almsB=read_fits_alms(fitsfile,3);
    disp('Loaded!')
elseif isstruct(fitsfile)
    almsT = fitsfile.almsT;
    almsE = fitsfile.almsE;
    almsB = fitsfile.almsB;
else
    error('Fits file should be either a struct or filename')
end


if isstr(cutparamfile)
    if length(cutparamfile)==4
        fname = ['data/real/' cutparamfile '_cutparams.mat'];
        if exist(fname,'file')
            cutparamfile = fname;
        else
            error('Year not found');
        end
    end
    
    disp('Loading cutparams file...')
    cutparams = load(cutparamfile);
    disp('Loaded!')
elseif isstruct(cutparamfile)
    cutparams = cutparamfile;    
else
    error('cutparamfile should be either a struct or filename')
end

if isstr(finalfile)
    disp('Loading final file...')
    final = load(finalfile) 
    disp('Loaded!')
elseif isstruct(finalfile)
    final = finalfile;
    clear finalfile
else
    error('Final file should be either a struct or filename')
end






%% Determine pair weights by coadding like we do with maps
% over phase, per pixel
% each tag will have per-pixel NET
disp('Calculating Pair Weights')

% find unique schedules
scheds = {};
tg = cutparams.cut_tags{1};
scheds{1} = tg(1,[1:8 12:17]);

for i = 2:length(cutparams.cut_tags)
    flag = true;
    tg = cutparams.cut_tags{i};
    tg = tg(1,[1:8 12:17]);
    for j = 1:length(scheds)
        if strcmp(tg,scheds{j})
            flag = false;
        end
    end
    if flag
        scheds{end+1} = tg;
    end
end

% build weights by phase
phases = {'B','C','E','F','G','H','I'};
nt = [];

tic;
for i = 1:length(scheds)
    for j = 1:length(phases)
        ind = [];
        for k = 1:10
            tname = scheds{i};
            tname = [tname(1:8) phases{j} sprintf('%02i',k) tname(9:end)];
            for m = 1:length(cutparams.cut_tags)
                if strcmp(tname,cutparams.cut_tags{m})
                    ind(end+1) = m;
                end
            end
        end
        
        nt(end+1,:) = wmean(cutparams.cpa.net(ind,:),cutparams.c2a.overall(ind,:),1);
    end
end
toc;

% average weights over all phases
ntperpair = mean(nt,1);
ntperpair(:,isnan(phiQ)) = NaN;
%% Rotate the alms by the phi angles of each channel.
disp('Rotating Spectra')
tic;

almsTd = almsT;
almsEd = almsE;
almsBd = almsB;

wtot = nansum(1./ntperpair);
almsEd.alms = [0 0];
almsBd.alms = [0 0];
for i = 1:length(phiQ)
    c = cosd(2*phiQ(i))./ntperpair(i)./wtot;
    s = sind(2*phiQ(i))./ntperpair(i)./wtot;

    almsEd.alms = almsEd.alms + cm*almsEd.alms+sm*almsBd.alms;
    almsBd.alms = almsBd.alms + cm*almsBd.alms+sm*almsEd.alms;
end

disp('.')
[l,TT]=alm2cl(almsTd);          disp('.')
[l,TE]=alm2cl(almsTd,almsEd);   disp('.')
[l,EE]=alm2cl(almsEd);          disp('.')
[l,BB]=alm2cl(almsBd);          disp('.')
[l,TB]=alm2cl(almsTd,almsBd);   disp('.')
[l,EB]=alm2cl(almsEd,almsBd);   disp('.')
[l,ET]=alm2cl(almsEd,almsTd);   disp('.')
[l,BT]=alm2cl(almsBd,almsTd);   disp('.')
[l,BE]=alm2cl(almsBd,almsEd);   disp('.')

%TT = 1; TE = 2; EE = 3;
%BB = 4; TB = 5; EB = 6;
%ET = 7; BT = 8; BE = 9;

%try to calc expectation values:
rotmod.l = l;
rotmod.C_l = abs([TT TE EE BB TB EB ET BT BE]);
lrep = repmat(l,1,9);
rotmod.Cs_l = lrep.*(lrep+1)./2./pi.*rotmod.C_l;


rrot = final.r;
rrot = calc_expvals(rrot, rotmod, final.bpwf);

finaldummy = final;
finaldummy.r.real = rrot.expv;
polopt.simtype = 'simd';
polopt.spec = {'EB'};


toc;
disp('fitting for rotation angle')
polrot = reduc_global_rotation(finaldummy, polopt, true,{'B2016comb'});


bins = 2:15;
freeparam.free = 1;
freeparam.lb = -90;
freeparam.ub = 90;
[alpha, alphaerr, gof, stat, cov] = matmin('chisq',...
            0, true,'EB_rotation_model',rrot.expv(bins,6), ones(length(bins),1)/0.1, rrot.expv(bins,3),rrot.expv(bins,4));



keyboard()

