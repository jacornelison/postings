function [mirr_model mirr_err] = rps_get_mirror_model(model,md,rpsopt,p)

fprintf('model params: [%2.3f %i]\n',model(1),model(2))


%%%%%%%%%%%%%%%%%%%%%%%%
% Fitting stuff
%Rearrange params
mount = rpsopt.mount;
%mount.aperture_offz = model(1);

mirror = rpsopt.mirror;
mirror.tilt = model(1);
mirror.roll = model(2);
source = rpsopt.source;

[mx, my, ex, ey] = deal([]);


prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;

bs.r = 0;
bs.theta = 0;
bs.chi = 0;
bs.chi_thetaref = 0;
bs.drumangle = 0;
mount.aperture_offr = 0; % Locate boresight 'aperture' in center of drum.



for i = 1:length(md.ch)
	
[r, theta, psi] = ...
    keck_beam_map_pointing(md.az{i}, md.el{i}, md.dk{i}, mount, ...
                           mirror, source, bs);

x = 2. * sind(r/2).*cosd(theta)*180.0/pi;
y = 2. * sind(r/2).*sind(theta)*180.0/pi;

% Trying matmin at Victor's suggestion
    % Found that it's infinitely better than fminsearch.
    % Build free parameter structure
    % [x0 y0 sigx sigy rho B0 psi xpol N1 N2 A]
    freepar.free = [1 1 1 1 1 1];
    freepar.lb = [-1e4 -1e3 -1e3 0 0 -1];
    freepar.ub = [1e4 1e3 1e3 2 2 1];
    
    % Concat timestream data together
    %[ang_data,x,y] = rps_concat_data(tods,ch);
    %x = tods{i}.x;
    %y = tods{i}.y;
    data = md.todcos{i};
    ind = abs(x-prx(md.ch(i))) < 100 & abs(y-pry(md.ch(i))) < 100;
    xi = x(ind);
    yi = y(ind);
    data = data(ind);
	[A_guess, mind] = max(md.todcos{i});
	guess = [A_guess prx(md.ch(i)) pry(md.ch(i)) p.fwhm_maj(md.ch(i))^2 p.fwhm_maj(md.ch(i))^2 0];
    % Very rough estimate of data rms noise by masking all data within 10sigma
    % of the beam, then calculate std.
    %dist = sqrt((x-guess(2)).^2+(y-guess(3)).^2);
    %ind = find(dist>mean(guess(4:5))*10);
    
    %keyboard()
    %sd = std(data(ind));
    %mn = mean(ang_data(ind));
    % Create model timestream
    
    
    %fprintf('Running sim %i of %i...\n',i,length(tods))
    
    % Estimate parameters
    %[param, err, gof, stat, cov] = matmin('chisq',guess, freepar,'egauss2', data, sd, x, y);
    
	%fun = @(z)chimin(z,data,x,y);
    fun = @(z)sum((data-egauss2(z,xi,yi)).^2);   
	[T, param] = evalc('fminsearch(fun,guess)');

    
	mx(end+1) = param(2);
	my(end+1) = param(3);
	
	%ex(end+1) = err(2);
	%ey(end+1) = err(3);
    %bparam.ch(i) = ch;
    %bparam.param(i,:) = param;
    %bparam.err(i,:) = err;
    %bparam.gof(i) = gof;
    %bparam.pte(i) = chi2cdf(gof,length(data)-length(guess));
    %bparam.stat(i) = stat;
    %bparam.dk(i) = tods{1}.scan.dk_ofs;
    %bparam.data_rms(i) = sd;
    
end

mirr_model = [mx'; my']; % Rotate to match prx/pry
%mirr_err = [ex'; ey'];

function chi2 = chimin(model,data,x,y)

chi2 = sum((data-egauss2(model,x,y)).^2);