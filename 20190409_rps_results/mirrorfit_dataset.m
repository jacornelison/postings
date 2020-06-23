
load('b3rpsfiles')
dirname = 'rps_data/2019_mirror/';

if strcmp(dirname(15:end-1),'azel')
    figfold = 'figs/';
elseif strcmp(dirname(15:end-1),'lockin')
    figfold = 'figs2/';
elseif strcmp(dirname(15:end-1),'new_decon')
    figfold = 'figs3/';
end


load('man_cuts')
good_chans = find(man_cuts);


%p = get_array_info(20180131,'ideal','ideal','ideal','ideal');
prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;
fd = [];
fd.bparam = [];
fd.berr = [];
fd.aparam = [];
fd.aerr = [];
fd.amodel = [];
fd.ch = [];
fd.sch = [];
fd.row = [];
fd.stat = [];
fd.dk = [];
fd.gof = [];
fd.rot = [];
fd.data_rms = [];
fd.inda = [];
fd.phi_s = [];
fd.phi_d = [];
fd.bchi2 = [];

% Manually cut tiles so far:
% [1 2 3 4 8 13 11]
rot = -180:30:180;
for i = 1:length(sch)
    schname = sch{i}.name(end-15:end-4);
    for j = 1:sch{i}.nrows
        rpsparam_file = [dirname, 'params/param_', schname, sprintf('_%02i.mat', j)];
        try
            load(rpsparam_file)
            if ~isempty(i_param) & ~isfield(i_param{1},'map_param')
                for k = 1:length(i_param)
                    if any(p_ind.rgl100==i_param{k}.ch)
                        if length(i_param{k}.rot)~=13
                            [bparam, berr] = deal(NaN(1,18));
                            bparam(1:5) = i_param{k}.bparam(1:5);
                            berr(1:5) = i_param{k}.berr(1:5);
                            drot = NaN(1,13);
                            ind = find(ismember(rot,round(i_param{k}.rot)));
                            
                            drot(ind) = i_param{k}.rot;
                            bparam(ind+5) = i_param{k}.bparam(6:end);
                            berr(ind+5) = i_param{k}.berr(6:end);
                            
                            i_param{k}.bparam = bparam;
                            i_param{k}.berr = berr;
                            i_param{k}.rot = drot;
                        end
                        
                        
                        fd.rot = [fd.rot; i_param{k}.rot];
                        fd.bparam = [fd.bparam; i_param{k}.bparam];
                        fd.berr = [fd.berr; i_param{k}.berr];
                        fd.amodel = [fd.amodel; rps_get_mod_model(...
                            i_param{k}.aparam, i_param{k}.rot)];
                        %                             if p.mce(i_param{k}.ch) == 0
                        %                                 i_param{k}.aparam(1) = i_param{k}.aparam(1) - 90;
                        %                                 i_param{k}.phi_d = i_param{k}.phi_d - 90;
                        %                                 i_param{k}.phi_s = i_param{k}.phi_s - 90;
                        %                             end
                        fd.aparam = [fd.aparam; i_param{k}.aparam];
                        fd.aerr = [fd.aerr; i_param{k}.aerr];
                        fd.ch = [fd.ch; i_param{k}.ch];
                        fd.sch = [fd.sch; i];
                        fd.row = [fd.row; j];
                        fd.stat = [fd.stat; i_param{k}.astat];
                        fd.dk = [fd.dk; i_param{k}.dk];
                        fd.gof = [fd.gof; i_param{k}.agof];
                        fd.phi_s = [fd.phi_s; i_param{k}.phi_s];
                        fd.phi_d = [fd.phi_d; i_param{k}.phi_d];
                        fd.bchi2 = [fd.bchi2; i_param{k}.bchi2];
                        %fd.gof = [fd.gof; sum(((fd.bparam(end,6:end)-fd.amodel(end,:))/i_param{k}.var).^2)];
                        
                        fd.data_rms = [fd.data_rms; i_param{k}.data_rms];
                        if strcmp(p.pol(i_param{k}.ch),'A')
                            if p.mce(i_param{k}.ch)~=0
                                fd.inda = [fd.inda; 1];
                            else
                                fd.inda = [fd.inda; 0];
                            end
                        else
                            if p.mce(i_param{k}.ch)~=0
                                fd.inda = [fd.inda; 0];
                            else
                                fd.inda = [fd.inda; 1];
                            end
                        end
                        
                    end
                end
            end
        catch exception
            s = 'couldNotReadFile';
            if strcmp(exception.identifier((end-length(s)+1):end),s)
                disp(['Could not read: ', rpsparam_file])
            else
                disp(exception.message)
                keyboard()
                disp(['Skipping: ', rpsparam_file])
            end
            continue
        end
    end
end

fd.norm_ang = atand(tand(fd.aparam(:,1)));
m = repmat(max(fd.bparam(:,6:end)'),[13,1])';
fd.res_var = var(fd.bparam(:,6:end)-fd.amodel(:,:),0,2);
fd.res_perc = (fd.bparam(:,6:end)-fd.amodel(:,:))./m;


% First, find some good detectors, start with one pair per focal plane
% We already have the tod's for them. We can do the quick and dirty way
% of manually selecting channels from the appropriate RPS TOD's.
% We'll want to do this for each dk as well.
% Perhaps start with a few channels to see how it does?
% Also, we don't need every single pol rotation, so that should save us some time.


dks = unique(fd.dk);

% row and column that has "good" data for each deck
col = [4 4 4 4 4 5 3 4 3 6 5 3 5 4 4 5 4 3 3 5]; %
row = [4 4 4 6 6 3 4 4 3 5 5 3 5 3 3 3 4 3 4 6]; %
%     [1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0]
pola =[4 4 4 4 4 4 4 1 4 4 4 1 1 4 4 4 1 1 4 4];
polb =[1 1 1 1 1 1 1 4 1 1 1 4 4 1 1 1 4 4 1 1];

md = []; %mirror data
md.todcos = {};
md.az     = {};
md.el     = {};
md.dk     = {};
md.dk0    = [];
md.ch     = [];

for j = 1:4
    % Choose which of the 13 TODs from the rasterset to load
    switch j
        case 1
            pola =[4 4 4 4 4 4 4 1 4 4 4 1 1 4 4 4 1 1 4 4];
            polb =[1 1 1 1 1 1 1 4 1 1 1 4 4 1 1 1 4 4 1 1];
        case 2
            pola =[4 4 4 4 4 4 4 1 4 4 4 1 1 4 4 4 1 1 4 4]+1;
            polb =[1 1 1 1 1 1 1 4 1 1 1 4 4 1 1 1 4 4 1 1]+1;
        case 3
            polb =[4 4 4 4 4 4 4 1 4 4 4 1 1 4 4 4 1 1 4 4];
            pola =[1 1 1 1 1 1 1 4 1 1 1 4 4 1 1 1 4 4 1 1];
        case 4
            polb =[4 4 4 4 4 4 4 1 4 4 4 1 1 4 4 4 1 1 4 4]+1;
            pola =[1 1 1 1 1 1 1 4 1 1 1 4 4 1 1 1 4 4 1 1]+1;
    end
    for i = [5 7 11 15 17]%1:20
        ch = find(p.tile == i & p.det_col==col(i) & p.det_row == row(i));
        %fprintf('tile %i: %i\n',i,length(fd.sch(find(fd.ch == ch(1)))))
        
        ind = find(fd.dk == dks(j) & fd.ch == ch(1));
        if length(ind)>1
            ind = ind(2); % Just in case there's 2 measurements at a single deck.
        end
        %fprintf('tile\t%i\tsch\t%i\trow\t%i\n',i,fd.sch(ind),fd.row(ind))
        if ~isempty(ind)
            schname = [sch{fd.sch(ind)}.name(end-15:end-4)];
            
            for ii = 1:2
                if ii == 1
                    fname = [dirname, 'tods/tod_', schname, sprintf('_%03i.mat', pola(i)+(fd.row(ind)-1)*13)];
                else
                    fname = [dirname, 'tods/tod_', schname, sprintf('_%03i.mat', polb(i)+(fd.row(ind)-1)*13)];
                end
                disp(['Loading ', fname])
                try
                load(fname)
                tod = rpstod.todcos(:,rpstod.ch==ch(ii));
                if size(tod,2)==1
                    md.todcos{end+1} = tod;
                    md.az{end+1} = rpstod.az;
                    md.el{end+1} = rpstod.el;
                    md.dk{end+1} = rpstod.dk;
                    md.dk0(end+1) = dks(j);
                    md.ch(end+1) = ch(ii);
                end
                catch
                    continue
                end
            end
        end
    end
end


%%%%%%%
% Setup fitting stuff


prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;

d = [prx(md.ch); pry(md.ch)];


% This shouldn't change. Ever.
rpsopt.source.azimuth = -177.6500060189; % +/- 4.5e-10
rpsopt.source.el = 2.3129580213; % +/- 1.6e-10
rpsopt.source.distance = 195.5; % +/- 1
rpsopt.source.height = rpsopt.source.distance*tand(rpsopt.source.el);

bs.r = 0;
bs.theta = 0;
bs.chi = 0;
bs.chi_thetaref = 0;
bs.drumangle = 0;
mount.aperture_offr = 0; % Locate boresight 'aperture' in center of drum.

rpsopt.mirror.tilt = 44.687;
rpsopt.mirror.roll = 0.148;

mm = rps_get_mirror_model([44.687 0.148], md, rpsopt, p);

tic
chisq([44.687 0.148],'rps_get_mirror_model',d,1,md,rpsopt,p)
toc

chi2 = [];
for i = 46.2:0.1:47.0
    mm = rps_get_mirror_model([i 0], md, rpsopt, p);
    chi2(end+1) = nansum(((d-mm)).^2);
end

fp.free = [1 1];
fp.ub = [46 0.2];
fp.lb = [44 -0.2];

guess = [44.687 0.148]; %Current best

[param, err, gof] = matmin('chisq',guess, fp,'rps_get_mirror_model',d,ones(size(d)),md,rpsopt,p);

%%

[mx,my] = deal([]);
fig = figure(1);
hold off
for i = 1:length(md.ch)
    %if ~any(ismember(md.ch(i),[705 710]))
    	
[r, theta, psi] = ...
    keck_beam_map_pointing(md.az{i}, md.el{i}, md.dk{i}, rpsopt.mount, ...
                           rpsopt.mirror, rpsopt.source, bs);

x = 2. * sind(r/2).*cosd(theta)*180.0/pi;
y = 2. * sind(r/2).*sind(theta)*180.0/pi;

    data = md.todcos{i};
    ind = abs(x-prx(md.ch(i))) < 100 & abs(y-pry(md.ch(i))) < 100;
    xi = x(ind);
    yi = y(ind);
    data = data(ind);
    [m,mind] = max(data);
    guess = [m prx(md.ch(i)) pry(md.ch(i)) p.fwhm_maj(md.ch(i))^2 p.fwhm_maj(md.ch(i))^2 0 mean(data)];
    fun = @(z)sum((data-egauss2(z,xi,yi)).^2);    
    %fun = @(z)chimin(z,tod,x,y);
    
    [T, param] = evalc('fminsearch(fun,guess)');
    
    mx(end+1) = param(2);
    my(end+1) = param(3);
    
    %mx(end+1) = xi(mind);
    %my(end+1) = yi(mind);
    
    if 1 
    hold off
    plot3(xi,yi,data)
    hold on
    plot3([1 1]*prx(md.ch(i)),[1 1]*pry(md.ch(i)),[0,1000],'r')
    %plot3([1 1]*xi(mind),[1 1]*yi(mind),[0,3000],'g')
    plot3([1 1]*param(2),[1 1]*param(3),[0,1000],'g')
    
    title(sprintf('channel %d, tile %d, row %d, col %d',md.ch(i),...
        p.tile(md.ch(i)), p.det_row(md.ch(i)), p.det_col(md.ch(i))))
        %pause(0.5)
    zzz = input(sprintf('dx %d dy %d',prx(md.ch(i))-param(2),pry(md.ch(i))-param(3)));
    end

end

if 1 
fig = figure(2);
clf(fig)
hold off
%plot(x,y)
hold on
%plot(mirror_model(1:end/2),mirror_model((end/2+1):end),'kx')
plot(mx,my,'gx')
plot(prx(md.ch),pry(md.ch),'rx')
%plot(prx,pry,'rx')
%plot(mx,my,'gx')
xlim([-20 20])
ylim([-20 20])
xlabel('x (^o)')
ylabel('y (^o)')
grid on
end




%% Get the target mirror params close using a single tod

% This shouldn't change. Ever.
rpsopt.source.azimuth = -177.6500060189; % +/- 4.5e-10
rpsopt.source.el = 2.3129580213; % +/- 1.6e-10
rpsopt.source.distance = 195.5; % +/- 1
rpsopt.source.height = rpsopt.source.distance*tand(rpsopt.source.el);

prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;

bs.r = 0;
bs.theta = 0;
bs.chi = 0;
bs.chi_thetaref = 0;
bs.drumangle = 0;
mount.aperture_offr = 0; % Locate boresight 'aperture' in center of drum.


% This can change.
rpsopt.mirror.tilt = 44.7;
rpsopt.mirror.roll =  0.12;


[r, theta, psi] = ...
    keck_beam_map_pointing(rpstod.az, rpstod.el, rpstod.dk, rpsopt.mount, ...
    rpsopt.mirror, rpsopt.source, bs);


x = 2 * sind(r / 2) .* cosd(theta) * 180 / pi;
y = 2 * sind(r / 2) .* sind(theta) * 180 / pi;



[mx,my] = deal([]);
fig = figure(1);
hold off
for i = 1:length(rpstod.ch)
    %if ~any(ismember(rpstod.ch(i),[705 710]))
    data = rpstod.todcos(:,rpstod.ch==rpstod.ch(i));
    ind = abs(x-prx(rpstod.ch(i))) < 2 & abs(y-pry(rpstod.ch(i))) < 2;
    xi = x(ind);
    yi = y(ind);
    data = data(ind);
    [m,mind] = max(data);
    %guess = [m prx(rpstod.ch(i)) pry(rpstod.ch(i)) p.fwhm_maj(rpstod.ch(i))^2 p.fwhm_maj(rpstod.ch(i))^2 0 mean(data)];
    %fun = @(z)sum((data-egauss2(z,xi,yi)).^2);
    
    %[T, param] = evalc('fminsearch(fun,guess)');
    
    %mx(end+1) = param(2);
    %my(end+1) = param(3);
    
    mx(end+1) = xi(mind);
    my(end+1) = yi(mind);
    
    if 1
        hold off
        plot3(xi,yi,data)
        hold on
        plot3([1 1]*prx(rpstod.ch(i)),[1 1]*pry(rpstod.ch(i)),[0,3000],'r')
        plot3([1 1]*xi(mind),[1 1]*yi(mind),[0,3000],'g')
        %plot3([1 1]*param(2),[1 1]*param(3),[0,3000],'g')
        
        title(sprintf('channel %d, tile %d, row %d, col %d',rpstod.ch(i),...
            p.tile(rpstod.ch(i)), p.det_row(rpstod.ch(i)), p.det_col(rpstod.ch(i))))
        %pause(0.5)
        zzz = input(sprintf('dx %d dy %d',prx(rpstod.ch(i))-param(2),pry(rpstod.ch(i))-param(3)));
    end
    
end


if 1
    fig = figure(2);
    clf(fig)
    hold off
    plot(x,y)
    hold on
    %plot(mirror_model(1:end/2),mirror_model((end/2+1):end),'kx')
    %plot(mx,my,'gx')
    plot(prx(rpstod.ch),pry(rpstod.ch),'rx')
    %plot(prx,pry,'rx')
    plot(mx,my,'gx')
    xlim([-5 5])
    ylim([-5 5])
    xlabel('x (^o)')
    ylabel('y (^o)')
    grid on
end





