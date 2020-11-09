function ad = rps_get_angle_data(dirname,p,rpsopt,PLOT)

if ~exist('PLOT','var')
    PLOT = false;
end

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
ad.param = [];
ad.ch = [];
ad.model = [];
ad.data = [];
ad.aparam = [];
ad.agof = [];
ad.astat = [];
ad.frac_res = [];
ad.res_var = [];
ad.dk = [];
% Build directory search
if fitopt.perdk 
    srchname = '*_dk_*.mat';
    ad.dk = [];
else
    srchname = '.mat';
end

m = 0;
fields = {'aparam','aerr', 'agof', 'astat', 'acov'};
flagged = false;

fname0 = [dirname 'params/angparam_ch_*_' fitopt.type];
D = dir([fname0 srchname]);
    
for i=1:length(D)


        try
            
            load([dirname 'params/' D(i).name]);
            
            ad.param{end+1,1} = parm;
            ad.aparam(end+1,:) = parm.aparam(1:6);
            ad.agof(end+1,1) = parm.agof;
            ad.astat(end+1,1) = parm.astat;
            ad.model{end+1,1} = model;
            ad.data{end+1,1} = data;
            
            ad.ch = [ad.ch; str2num(D(i).name(13:16))];
            
            if fitopt.perdk
                dkind = 22+length(fitopt.type);
                ad.dk = [ad.dk; str2num(D(i).name(dkind:dkind+2))];
            end
            
            ad.frac_res{end+1,1} = (model-data)/max(data);
            ad.res_var(end+1,1) = nanvar(ad.frac_res{end});
            l = length(parm.aparam(7:end));
            if l>m
                m = l;
            end
        catch exception
            switch exception.identifier
                case 'MATLAB:dimagree'
                    if ~flagged
                        disp('rcvd dimension error. Skipping.')
                    end
                    flagged = true;
                    continue
                otherwise
                    disp(exception.message)
                    keyboard()
                    continue
             end
            
            
        end

end

