function [source_model, mx, my] = rps_get_source_model(model,md,moonopt,p)

fprintf('model params: [%5.5f %5.5f %5.5f]\n',model(1),model(2),model(3))


%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting stuff
%Rearrange params
mount = moonopt.mount;
%mount.aperture_offz = model(1);

mirror = moonopt.mirror;
%source = moonopt.source;

source.azimuth = model(1);
source.el = model(2);
source.distance = model(3);
source.height = source.distance*tand(source.el);


[mx, my, ex, ey] = deal([]);


prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;

bs.r = 0;
bs.theta = 0;
bs.chi = 0;
bs.chi_thetaref = 0;
bs.drumangle = 0;
mount.aperture_offr = 0; % Locate boresight 'aperture' in center of drum.

%keyboard()
% 
% [r, theta, psi] = ...
%     keck_beam_map_pointing(md.az{i}, md.el{i}, md.dk{i}, mount, ...
%                            mirror, source, bs);
%     x0 = 2. * sind(r/2).*cosd(theta)*180.0/pi;
%     y0 = 2. * sind(r/2).*sind(theta)*180.0/pi;
%                        
for i = 1:length(md.ch)
	
    % Don't want to have to recompute this every time. Just when the az,el
    % changes.
    if i == 1 | ((md.ind(i-1)-md.ind(i))<0)
        [r, theta, psi] = ...
            keck_beam_map_pointing(md.az{i}, md.el{i}, md.dk{i}, mount, ...
            mirror, source, bs);
        x0 = 2. * sind(r/2).*cosd(theta)*180.0/pi;
        y0 = 2. * sind(r/2).*sind(theta)*180.0/pi;
    end
    
    data = md.tod{i};

    % Lots of data, so isolate to just 3degrees around beam center
    ind = abs(x0-prx(md.ch(i))) < 5 & abs(y0-pry(md.ch(i))) < 5;
    
    x = x0(ind);
    y = y0(ind);
    data = data(ind);
    
    % Trying matmin at Victor's suggestion
    % Found that it's infinitely better than fminsearch.
    % Build free parameter structure
    freepar.free = [1 1 1 1 1 1];
    freepar.lb = [-1e4 -1e3 -1e3 0 0 -1];
    freepar.ub = [1e4 1e3 1e3 2 2 1];
    
    
    [A_guess, mind] = max(data);
    
    if strcmp(moonopt.fittype,'fminsearch')
        guess = [A_guess prx(md.ch(i)) pry(md.ch(i)) p.fwhm_maj(md.ch(i))^2 p.fwhm_maj(md.ch(i))^2 0];
        fun = @(z)chimin(z,data,x,y);
        [T, param] = evalc('fminsearch(fun,guess)');
        
        mx(end+1) = param(2);
        my(end+1) = param(3);
    elseif strcmp(moonopt.fittype,'max')
        % If the there's no data for a certain channel,
        % that means the mirror params are way off, so we should make the
        % offset super high to reflect a bad fit.
        if ~isempty(mind)
            mx(end+1) = x(mind);
            my(end+1) = y(mind);
        else
            mx(end+1) = 1e4;
            my(end+1) = 1e4;
        end
    elseif strcmp(moonopt.fittype,'matmin')
        % Estimate parameters
        freepar.free = [1 1 1 1 1 1];
        freepar.lb = [-1e4 -1e3 -1e3 0 0 -1];
        freepar.ub = [1e4 1e3 1e3 2 2 1];
    
        guess = [A_guess prx(md.ch(i)) pry(md.ch(i)) p.fwhm_maj(md.ch(i))^2 p.fwhm_maj(md.ch(i))^2 0];
        [param, err, gof, stat, cov] = matmin('chisq',guess, freepar,'egauss2', data, ones(size(data)), x, y);
        mx(end+1) = param(2);
        my(end+1) = param(3);
    end
    
    % Very rough estimate of data rms noise by masking all data within 10sigma
    % of the beam, then calculate std.
    %dist = sqrt((x-guess(2)).^2+(y-guess(3)).^2);
    %ind = find(dist>mean(guess(4:5))*10);
    %sd = std(data(ind));
    %mn = mean(ang_data(ind));
    % Create model timestream
    
    
    %fprintf('Running sim %i of %i...\n',i,length(tods))
    
    % Estimate parameters
    %[param, err, gof, stat, cov] = matmin('chisq',guess, freepar,'egauss2', data, sd, x, y);
end

source_model = [mx'; my']; % Rotate to match prx/pry
%mirr_err = [ex'; ey'];

function chi2 = chimin(model,data,x,y)

chi2 = sum((data-egauss2(model,x,y)).^2);