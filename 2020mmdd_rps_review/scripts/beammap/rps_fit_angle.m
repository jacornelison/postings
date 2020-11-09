function [parm model data] = rps_fit_angle(sch,p,rpsopt,dirname,fd,ch)
% function a_param = rps_fit_angle(sch,p,rpsopt,fd,ch)
% Per channel fit to mod curves across all deck angles to account for
% precession and mispointing.
fprintf('\nFitting channel %04i\n.\n.\n.',ch)

if ~isfield(rpsopt,'fitopt')
    fitopt = [];
else
    fitopt = rpsopt.fitopt;
end

if ~isfield(fitopt,'type')
    fitopt.type = 'var_amp';
end

% Fit all data across DK-angles by default.
if ~isfield(fitopt,'perdk')
   fitopt.perdk = false;
elseif fitopt.perdk & ~strcmp(fitopt.type,'custom')
   fitopt.type = 'const_amp';
end

if isempty(find(fd.ch==ch))
    fprintf('No fit data for channel %04i detected. Exiting...\n.\n.\n.',ch);
    return
end

nsch = fd.sch(fd.ch==ch);
unsch = unique(nsch);
nrow = fd.row(fd.ch==ch);
ndk = fd.dk(fd.ch==ch);

[A, B3, B1] = deal(NaN(length(nsch),3));
[data, berr] = deal(NaN(length(nsch),18-5));
amp_guess = NaN(1,length(nsch));
count  = 1;

% Collect relevant per DK information
for i = 1:length(unsch)
    % Init vars
    ind = find(fd.sch==unsch(i) & fd.ch==ch);
    
    % Start by loading the tods specific to that channel
    nrps = sch{nsch(i)}.nrps;
    tods = rps_load_tods(sch,unsch(i),[1:nrps]+nrps*(nrow(i)-1),rpsopt,p,dirname,ch);
    
    % Concat the data
    [tod,x,y,phi,az,el,dk, rot] = rps_concat_data(tods,ch);
    
    % Acquire boresight position, ort and pointing timestreams.
    [A0,B3_0,B1_0] = kbmp_mount(az,el,dk+unique(p.drumangle),rpsopt.mount,0);
    
    % Reflect with mirror
    [mpos, mnorm] = kbmp_mirror(az,el,rpsopt.mount,rpsopt.mirror);
    [A0,B3_0,B1_0] = kbmp_reflect(A0, B3_0, B1_0, mpos, mnorm);
    % If there's more than one measurement in a single schedule (rare)
    for j = 1:length(ind) 
        bparam = fd.bparam(ind(j),:);
        data(count,:) = fd.bparam(ind(j),6:end);
        berr(count,:) = fd.berr(ind(j),6:end);
        amp_guess(1,count) = nanmax(data(count,:));
        % Find points that correspond with estimated beam centers.
        for k = 1:3
            A(count,k) = griddata(x,y,A0(:,k),bparam(1),bparam(2));
            B3(count,k) = griddata(x,y,B3_0(:,k),bparam(1),bparam(2));
            B1(count,k) = griddata(x,y,B1_0(:,k),bparam(1),bparam(2));
        end
        count = count+1;
    end
end


%% Fit stuff
% Use ideal pol angle as initial guess
ang_guess = atand(tand(p.chi(ch)+p.chi_thetaref(ch)));

% RPS mispointing is in local RPS az/el.
% Assume RPS is pointing straight at telescope.
%point_guess = [0 -atand(rpsopt.source.height/rpsopt.source.distance)];

% Assume RPS is completely level WRT to local gravity. I.e. pointing along
% the positive elevation.
point_guess = [0 atand(rpsopt.source.height/rpsopt.source.distance)];
coll_guess = [0 0];


% Different parameter options. This may change the number of output params
switch fitopt.type
    case 'const_amp'
        amp_guess = nanmean(amp_guess);
        freepar.free = [1 1 1 1 1 1 ones(1,length(amp_guess))];
        freepar.lb = [-10+ang_guess -1 -15 -15 -15 -15 zeros(1,length(amp_guess))];
        freepar.ub = [10+ang_guess 1 15 15 15 15 zeros(1,length(amp_guess))+1e6];
        guess = [ang_guess 0 coll_guess point_guess amp_guess];
    case 'var_amp'
        guess = [ang_guess 0 coll_guess point_guess amp_guess];
        freepar.free = [1 1 1 1 1 1 ones(1,length(amp_guess))];
        freepar.lb = [-10+ang_guess -1 -15 -15 -15 -15 zeros(1,length(amp_guess))];
        freepar.ub = [10+ang_guess 1 15 15 15 15 zeros(1,length(amp_guess))+1e6];
        
    case 'custom'
        
        

        if isfield(fitopt,'guess')
            guess = fitopt.guess;
        end
        if isfield(fitopt,'point_guess')
            point_guess = fitopt.point_guess;
        end
        if isfield(fitopt,'coll_guess')
            coll_guess = fitopt.coll_guess;
        end
        if isfield(fitopt,'ang_guess')
            ang_guess = fitopt.ang_guess;
        end
        if isfield(fitopt,'amp_guess')
            amp_guess = fitopt.amp_guess;
        elseif fitopt.perdk
            amp_guess = nanmean(amp_guess);
        end
        
        guess = [ang_guess 0 coll_guess point_guess amp_guess];
        freepar.free = [1 1 1 1 1 1 ones(1,length(amp_guess))];
        freepar.lb = [-10+ang_guess -1 -15 -15 -15 -15 zeros(1,length(amp_guess))];
        freepar.ub = [10+ang_guess 1 15 15 15 15 zeros(1,length(amp_guess))+1e6];
        
        if isfield(fitopt,'free')
            freepar.free = fitopt.free;
        end
        if isfield(fitopt,'lb')
            freepar.lb = fitopt.lb;
        end
        if isfield(fitopt,'ub')
            freepar.ub = fitopt.ub;
        end
                
    otherwise
        error('fit-type unknown')
end

fname0 = [dirname sprintf('params/angparam_ch_%04i_',ch) fitopt.type];
data0 = data;
if fitopt.perdk
    for i = 1:length(ndk)
        data = data0(i,:);
        nndk = sum(ndk(1:i)==ndk(i));
        [parm.aparam, parm.aerr, parm.agof, parm.astat, parm.acov] = matmin('chisq',...
            guess, freepar,	'rps_get_mod_model_vectors',data,berr(i,:),rot,A(i,:),B3(i,:),B1(i,:),rpsopt.source);
        
        model = rps_get_mod_model_vectors(parm.aparam,rot,A(i,:),B3(i,:),B1(i,:),rpsopt.source);
        fname = [fname0 sprintf('_dk_%03i_%i',floor(ndk(i)-rem(ndk(i),1)),nndk)];
        save(fname,'parm','model','data')
        disp(['Saved file to:' fname]);
    end
else
    data = reshape(data,1,[]);
    berr = reshape(berr,1,[]);

    [parm.aparam, parm.aerr, parm.agof, parm.astat, parm.acov] = matmin('chisq',...
        guess, freepar,	'rps_get_mod_model_alldk',data,berr,rot,A,B3,B1,rpsopt.source);
    
    fname = fname0;
    model = rps_get_mod_model_alldk(parm.aparam,rot,A,B3,B1,rpsopt.source);
    save(fname,'parm','model','data')
    disp(['Saved file to:' fname]);
end


fprintf('Angle fit of channel %04i complete!\n.\n.\n.\n',ch)
