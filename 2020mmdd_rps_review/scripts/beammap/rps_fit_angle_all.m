function [parm, cdata, model] = rps_fit_angle_all(sch,p,rpsopt,dirname,fd,chans,id)
% function rps_fit_angle_all(sch,p,rpsopt,dirname,fd,chans)
% Fit over multiple channels which may have multiple mod curves
% fit params: [phi_d 1-xpol E H Az_s El_s G]
%   phi_d, 1-xpol are per-channel
%   G is per mod-curve?
%   E, H, Az_s, El_s - just one per fit
%
% Top-level: Collect relevant info. Set up for fitting.
% 2nd: break up input data to use on rps_get_mod_model_vectors per channel,
% per measurement.
% 3rd: this should be rps_get_mod_model_vectors
%

fprintf('\nFitting over %04i channels \n.\n.\n.\n',length(chans))

source = rpsopt.source;

if ~isfield(rpsopt,'fitopt')
    fitopt = [];
else
    fitopt = rpsopt.fitopt;
end

if ~isfield(fitopt,'type')
    fitopt.type = 'default';
end

% We need arrays of data (bparam), errors (berr), A, B3, and B1.
% We should also track the dk and channel each index corresponds to.
fields = {'data','err','A','B3','B1','ch','dk'};

for i = 1:length(fields)
    cdata.(fields{i}) = [];
end

amp_guess = [];
chans = sort(chans); % Just so that we know for sure what order they are.

% Build Channel Data Struct.
for h = 1:length(chans)
    
    ch = chans(h);
    if isempty(find(fd.ch==ch))
        fprintf('\nNo fit data for channel %04i detected. Skipping\n.\n.\n.',ch);
        continue
    end
    
    nsch = fd.sch(fd.ch==ch);
    unsch = unique(nsch);
    nrow = fd.row(fd.ch==ch);
    ndk = fd.dk(fd.ch==ch);

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
            
            len = length(cdata.ch)+1;
            cdata.A(len,1:3) = NaN(1,3);
            cdata.B3(len,1:3) = NaN(1,3);
            cdata.B1(len,1:3) = NaN(1,3);
            
            % Find points that correspond with estimated beam centers.
            for k = 1:3
                %A(count,k) = griddata(x,y,A0(:,k),bparam(1),bparam(2));
                %B3(count,k) = griddata(x,y,B3_0(:,k),bparam(1),bparam(2));
                %B1(count,k) = griddata(x,y,B1_0(:,k),bparam(1),bparam(2));
                
                cdata.A(len,k) = griddata(x,y,A0(:,k),bparam(1),bparam(2));
                cdata.B3(len,k) = griddata(x,y,B3_0(:,k),bparam(1),bparam(2));
                cdata.B1(len,k) = griddata(x,y,B1_0(:,k),bparam(1),bparam(2));
            end
            
            cdata.data(len,:) = fd.bparam(ind(j),6:end);
            cdata.err(len,:) = fd.berr(ind(j),6:end);
            cdata.ch(len) = ch;
            cdata.dk(len) = fd.dk(ind(j));
        end
    end
end


%% Fit stuff
% rps_get_mod_model_vectors(aparam,eta,A,B3,B1,source,MIRROR,PLOT)
% 

% Angle Guesses - one per channel
phi_guess = p.chi(chans)+p.chi_thetaref(chans);
phi_guess = atand(tand(phi_guess'-45))+45; % put angles around 0 and 90

% Xpol efficiency - one per channel
xpol_guess = zeros(1,length(chans));

% Nutation angles - one for whole fit
if isfield(fitopt,'nut_guess')
nut_guess = fitopt.nut_guess; % [E H]
else
nut_guess = [0 0]; % [E H]
end
% misalignment angles - one for whole fit
if isfield(fitopt,'aln_guess')
aln_guess = fitopt.aln_guess; % [E H]
else
aln_guess = [0 0]; % [E H]
end

% Gain - one per measurements
g_guess = nanmax(cdata.data,[],2)';

p1 = ones(1,length(phi_guess));
x1 = ones(1,length(xpol_guess));
g1 = ones(1,length(g_guess));

if strcmp(fitopt.type,'fixed')
    n1 = zeros(1,length(nut_guess));
    a1 = zeros(1,length(aln_guess));
else
    n1 = ones(1,length(nut_guess));
    a1 = ones(1,length(aln_guess));
end

freepar.free = [p1 x1 n1 a1 g1];
freepar.lb = [-15+phi_guess -1*x1 -45*n1 -45*a1 0*g1];
freepar.ub = [15+phi_guess x1 45*n1 45*a1 1e4*g1];
guess = [phi_guess xpol_guess nut_guess aln_guess g_guess];

data = reshape(cdata.data,1,[]);
berr = reshape(cdata.err,1,[]);
rot = -180:30:180;

% For testing the model
if 0
    keyboard()
    tic; datatest = rps_get_mod_model_allch(guess,rot,cdata,source); toc;

    tic;
    [parm.aparam, parm.aerr, parm.agof, parm.astat, parm.acov] = matmin('chisq',...
        guess, freepar,	'rps_get_mod_model_allch',datatest,berr,rot,cdata,rpsopt.source);
    toc;

    tic; model = rps_get_mod_model_allch(parm.aparam,rot,cdata,source); toc;
end

tic;
[parm.aparam, parm.aerr, parm.agof, parm.astat, parm.acov] = matmin('chisq',...
    guess, freepar,	'rps_get_mod_model_allch',data,berr,rot,cdata,rpsopt.source);
toc;

model = rps_get_mod_model_allch(parm.aparam,rot,cdata,source);

if length(chans==1)
    fname = [dirname 'params/aparam_allfit_ch_' num2str(chans)];
else
    fname = [dirname 'params/aparam_allfit_' num2str(id)];
end

save(fname,'parm','cdata','guess','model')

disp(['Complete! Saving to:' fname])
pause(1)

