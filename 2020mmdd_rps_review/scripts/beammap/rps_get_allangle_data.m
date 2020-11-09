function ad = rps_get_allangle_data(dirname,p,rpsopt,PERCHAN)

if ~exist('PERCHAN','var')
    PERCHAN = false;
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
ad.phi = [];
ad.cdata = [];
ad.phi = [];
ad.xpol = [];
ad.nut = [];
ad.align = [];


flagged = false;
if PERCHAN
    fname0 = [dirname 'params/aparam_allfit_ch*'];
    D = dir(fname0);
else
    fname0 = [dirname 'params/aparam_allfit*'];
    D0 = dir(fname0);
    
    for i = 1:length(D0)
        if ~strcmp(D0(i).name(15:16),'ch')
            D(i) = D0(i);
        end
    end
    
end

    
for i=1:length(D)


        try
            
            load([dirname 'params/' D(i).name]);
            chans = unique(cdata.ch)';
            lp = length(chans);
            lx = lp;
            ln = 2;
            la = 2;
            lg = size(cdata.data,1);
            
            ad.ch = [ad.ch; chans];
            ad.phi = [ad.phi; parm.aparam(1:lp)'];
            ad.xpol = [ad.xpol; parm.aparam(lp+(1:lx))'];
            
            % ad.gain 
            
            for j = 1:length(chans)
            ad.nut(end+1,1:ln) = parm.aparam(lp+lx+(1:ln));
            ad.align(end+1,1:la) = parm.aparam(lp+lx+ln+(1:la));
            ad.param{end+1,1} = parm;
            ad.agof(end+1,1) = parm.agof;
            ad.astat(end+1,1) = parm.astat;
            ad.model{end+1,1} = model;
            ad.data{end+1,1} = reshape(cdata.data,1,[]);
            ad.cdata{end+1,1} = cdata;
            ad.frac_res{end+1,1} = (model-ad.data{end})/max(ad.data{end});
            ad.res_var(end+1,1) = nanvar(ad.frac_res{end});
            end
            
        catch exception
            switch exception.identifier
                case 'MATLAB:dimagree'
                    
                    if ~flagged
                        fprintf(['\nrcvd dimension error at line: ' num2str(exception.stack.line) '\n'])
                    end
                    keyboard()
                    flagged = true;
                    continue
                otherwise
                    fprintf(['\n' exception.message '\nat line: ' num2str(exception.stack.line) '\n'])
                    keyboard()
                    continue
             end
            
            
        end

end

