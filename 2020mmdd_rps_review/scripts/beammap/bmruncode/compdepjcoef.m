function compdepjcoef(compopt)
% function compdepjcoef(compopt)
%
% Compare deprojection coefficients (or beam map-derived values) 
% from any combination of:
% real data maps/beam map sims/fits from beam maps
% 
% INPUTS
%   
%   expt:       'bicep2','keck'
%   year:       2012, 2013, 2014, 2015
%   tocompare:  cell array containing the two things to compare:
%               if a map, loads it and scales dp coefficients
%               e.g. {'maps/1450/real_a_filtp3_weight3_gs_dp1111', ...
%                     'beammap'} uses values from 'beamwid' file
%               The first is on the x axis, second on the y
%   lab:        cell array indicating the two things
%               to compare - must be one of the 3 following:
%               'bm_meas','bm_dp','real_dp'
%   perwhat:    for Keck, do you want the plots
%               peryear: all combined
%               perrx:   each rx individually
%               perfreq: combine over freqs

% Parse compopt
expt = compopt.expt;
year = compopt.year;
tocompare = compopt.tocompare;
if ~isfield(compopt,'lab')
  lab = {'bm_meas','real_dp'}; 
else
  lab = compopt.lab;
end
if ~isfield(compopt,'perwhat')
  perwhat = 'peryear';
else
  perwhat = compopt.perwhat;
end

chflags = get_default_chflags(expt,num2str(year));

% Get ukpv (for scaling) and beam map-derived values
switch expt
  case 'bicep2'
    [p,ind] = get_array_info('20120505','obs','obs','obs','obs',...
                             chflags,'obs','obs');
    p = rmfield(p,'expt');
    [p1,kk] = ParameterRead(['aux_data/beams/' ...
                             'beams_bicep2_obs_rwa_20130607.csv']);
  case 'keck'
    [p,ind] = get_array_info([num2str(year) '0505'],'obs','obs',...
                             'obs','obs',chflags,'obs','obs');
    p = rmfield(p,'expt');
    [p1,kk] = ParameterRead(['aux_data/beams/' ...
                             'beamwid_' num2str(year) '0101.csv']);
end

% Following based on /n/home06/csheehy/papers/bicep2_systematics/code

% Get coefficients into mutually compatible formats - keep them in full
% length (i.e. per-detector) for now and then downselect later
for ii = 1:2
  switch tocompare{ii}
    case 'beammap'
      coeff(ii).b = massage_bmdata(p,ind,p1);
    otherwise
      m = load(tocompare{ii},'coaddopt');
      coeff(ii).b = massage_dpcoeff(m.coaddopt.b,m.coaddopt.bw,...
                                    year,p,ind); 
  end
end


% Figure out how to divide things up and plot
switch perwhat
  case 'peryear'
    make_dpplot(coeff(1).b,coeff(2).b,lab, ...
                expt,year,perwhat,[]);    
  case 'perrx'
    rxs = unique(p.rx);
    for ii = 1:length(rxs)
      idx = find(p.rx == rxs(ii));
      make_dpplot(coeff(1).b(idx,:),coeff(2).b(idx,:), ...
                  lab,expt,year,perwhat,['Rx' num2str(rxs(ii)) '_']);
    end
  case 'perfreq'
    freqs = unique(p.band(p.band ~= 0));
    for ii = 1:length(freqs)
      idx = find(p.band == freqs(ii));
      make_dpplot(coeff(1).b(idx,:),coeff(2).b(idx,:), ...
                  lab,expt,year,perwhat,[num2str(freqs(ii)) '_']);
    end
end


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coeffout = massage_dpcoeff(b,bw,year,p,ind)
% Take the 'b' matrix from a coadded map and make it, like, useful - 
% weighted average over time
% b = map.coaddopt.b, bw = map.coaddopt.bw
% where b = size(numel(ind.a),2,numel(dp coeffs),....)
% coeffout should be ndets x n_dp coeffs, calibrated 
% We'll cut it down later

% Expand weights in bw to size of b
bw = reshape(bw,[size(bw,1),size(bw,2),1,size(b,4)]);
bw = repmat(bw,[1,1,6,1]); % 6 = number of dp coeffs - could be more?

% Average over left-right scans
b = squeeze(nansum(b.*bw,2)./nansum(bw,2));
bw = squeeze(nansum(bw,2));

% Average over time
b = nansum(b.*bw,3)./nansum(bw,3);

% Calibrate
% We index 'b' with a-based indices (e.g. length 256, 1056, 1280)
% So find the indices OF ind.a which are RGLs
[dum rglind] = intersect(ind.a,ind.rgla);
% Get sigma and ukpv, which scale the coefficients, to the same length
sigma = ((p.fwhm_maj + p.fwhm_min)/2)/(2*sqrt(2*log(2)));
sigma = sigma(ind.a);
ukpv = p.ukpv(ind.a);

% Do this per-pair - why not?
b(rglind,1) = b(rglind,1)*2.*ukpv(rglind);
b(rglind,2) = b(rglind,2)*2.*ukpv(rglind);
b(rglind,3) = b(rglind,3)*2.*ukpv(rglind);
b(rglind,4) = b(rglind,4)*2.*ukpv(rglind)./sigma(rglind);
b(rglind,5) = b(rglind,5).*4.*ukpv(rglind)./(sigma(rglind).^2);
b(rglind,6) = b(rglind,6).*4.*ukpv(rglind)./(sigma(rglind).^2);

% NaN out anything that's identically zero
b(b == 0) = NaN;

% This is currently in n_pairs, so bring it up to n_dets
coeffout = NaN(length(ind.e),6);
coeffout(ind.a,1) = b(:,1);
coeffout(ind.a,2) = b(:,2);
coeffout(ind.a,3) = b(:,3);
coeffout(ind.a,4) = b(:,4);
coeffout(ind.a,5) = b(:,5);
coeffout(ind.a,6) = b(:,6);

% Rotate to account for drum angle
% dx/dy
[theta rho] = cart2pol(coeffout(:,2),coeffout(:,3));
[coeffout(:,2) coeffout(:,3)] = ...
    pol2cart(theta + p.drumangle(ind.e)*pi/180,rho);
% dp/dc
R(1,1,:) = cosd(2*p.drumangle(ind.e));
R(1,2,:) = -sind(2*p.drumangle(ind.e));
R(2,1,:) = sind(2*p.drumangle(ind.e));
R(2,2,:) = cosd(2*p.drumangle(ind.e));
tmp = [coeffout(:,5)'; coeffout(:,6)'];
tmpnew = NaN(size(tmp));
for ii = 1:size(tmp,2)
  tmpnew(:,ii) = R(:,:,ii)*tmp(:,ii);
end
coeffout(:,5) = tmpnew(1,:);
coeffout(:,6) = tmpnew(2,:);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function coeffout = massage_bmdata(pp,ind,p1)
% Take parameters from get_array_info and beamcen files 
% (in principle with enough 'obs' tags should be the same)
% and transform into values comparable to dp coefficients
% pp from get_array_info - for relgain
% p1 from beamcen - for everything else
% coeffout should be npairs x n_dp coeffs, calibrated

coeffout = NaN(length(ind.e),6); 

% Take gain from get_array_info values
g = 1./pp.ukpv;
gbar = (g(ind.rgla) + g(ind.rglb))/2;
abg = (g(ind.rgla) - g(ind.rglb))./gbar;

% Other beam params from get_array_info - if using this would need to
% account for drum angle too
%{
[s c p] = egauss2_mmt2scp(pp.fwhm_maj,pp.fwhm_min,pp.alpha + pp.theta);
% Diff pointing for horizon dk = 0
[dec ra] = reckon(0,0,pp.r,pp.theta - 90 - 180);
% Use -dec here, so +y -> az = -dec
[x y] = radec_to_gnomonic(ra,-dec,0,0);
dx = x(ind.rgla) - x(ind.rglb);
dy = y(ind.rgla) - y(ind.rglb);
% Other diff params
ds = s(ind.rgla) - s(ind.rglb);
dc = c(ind.rgla) - c(ind.rglb);
dp = p(ind.rgla) - p(ind.rglb);
%}

% From beamcen file, take differences for dx/dy/dsig
p1.x = p1.r.*cosd(p1.theta);
p1.y = p1.r.*sind(p1.theta);
dx = p1.x(ind.rgla) - p1.x(ind.rglb);
dy = p1.y(ind.rgla) - p1.y(ind.rglb);
ds = p1.sigma(ind.rgla) - p1.sigma(ind.rglb);

% If pre-calculated diff ellip exist, use them:
if isfield(p1,'dp')
  dp = p1.dp(ind.rgla);
  dc = p1.dc(ind.rgla);
else % Otherwise take the differences here
  dp = p1.p(ind.rgla) - p1.p(ind.rglb);
  dc = p1.c(ind.rgla) - p1.c(ind.rglb);
end

coeffout(ind.rgla,1) = abg;
coeffout(ind.rgla,2) = dx;
coeffout(ind.rgla,3) = dy;
coeffout(ind.rgla,4) = ds;
coeffout(ind.rgla,5) = dp;
coeffout(ind.rgla,6) = dc;

% NaN out anything that's identically zero
coeffout(coeffout == 0) = NaN;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_dpplot(coeff1,coeff2,lab,expt,year,perwhat,rxfreq)
% 6-panel standard plot - make both density and scatter plots
% Need to add bias line

figure(1); clf;
setwinsize(gcf,400,600)

figure(2); clf;
setwinsize(gcf,400,600)

nb = 20; % 20 bins seems nice

% Relgain
figure(1)
subplot(3,2,1)
x = cvec(coeff1(:,1));
y = cvec(coeff2(:,1));
cutnan = find(isnan(x)); x(cutnan) = []; y(cutnan) = [];
cutnan = find(isnan(y)); x(cutnan) = []; y(cutnan) = [];
P = corrcoef(x,y);
PP = polyfit(x,y,1);
lx = -0.05;
hx = 0.05;
[xt yt n] = hfill2(x,y,nb,lx,hx,nb,lx,hx);
imagesc(xt,yt,n);
set(gca,'YDir','normal');
axis square
xlim([lx hx]);
ylim(xlim());
hold on;
plot([-1 1],[-1 1],'k');
text(lx+(hx-lx)*.05,hx-(hx-lx)*.05,'\deltag',...
     'HorizontalAlignment','Left','VerticalAlignment','top',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.15,...
     ['r=' sprintf('%0.2f',P(1,2))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.05,...
     ['m=' sprintf('%1.2f',PP(1))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);

figure(2)
subplot(3,2,1)
scatter(coeff1(:,1),coeff2(:,1),'.');
axis square
xlim([lx hx]);
ylim(xlim());
hold on;
plot([-1 1],[-1 1],'k');
text(lx+(hx-lx)*.05,hx-(hx-lx)*.05,'\deltag',...
     'HorizontalAlignment','Left','VerticalAlignment','top',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.15,...
     ['r=' sprintf('%0.2f',P(1,2))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.05,...
     ['m=' sprintf('%1.2f',PP(1))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
clear x y P PP

% dx
figure(1)
subplot(3,2,3)
x = cvec(coeff1(:,2));
y = cvec(coeff2(:,2));
cutnan = find(isnan(x)); x(cutnan) = []; y(cutnan) = [];
cutnan = find(isnan(y)); x(cutnan) = []; y(cutnan) = [];
P = corrcoef(x,y);
PP = polyfit(x,y,1);
lx = -0.04;
hx = 0.04;
[xt yt n] = hfill2(x,y,nb,lx,hx,nb,lx,hx);
imagesc(xt,yt,n);
set(gca,'YDir','normal');
axis square
xlim([lx hx]);
ylim(xlim());
hold on;
plot([-1 1],[-1 1],'k');
text(lx+(hx-lx)*.05,hx-(hx-lx)*.05,'\deltax (deg)',...
     'HorizontalAlignment','Left','VerticalAlignment','top',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.15,...
     ['r=' sprintf('%0.2f',P(1,2))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.05,...
     ['m=' sprintf('%1.2f',PP(1))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);

figure(2)
subplot(3,2,3)
scatter(coeff1(:,2),coeff2(:,2),'.');
axis square
xlim([lx hx]);
ylim(xlim());
hold on;
plot([-1 1],[-1 1],'k');
text(lx+(hx-lx)*.05,hx-(hx-lx)*.05,'\deltax (deg)',...
     'HorizontalAlignment','Left','VerticalAlignment','top',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.15,...
     ['r=' sprintf('%0.2f',P(1,2))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.05,...
     ['m=' sprintf('%1.2f',PP(1))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
clear x y P PP

% dy
figure(1)
subplot(3,2,5)
x = cvec(coeff1(:,3));
y = cvec(coeff2(:,3));
cutnan = find(isnan(x)); x(cutnan) = []; y(cutnan) = [];
cutnan = find(isnan(y)); x(cutnan) = []; y(cutnan) = [];
P = corrcoef(x,y);
PP = polyfit(x,y,1);
lx = -0.04;
hx = 0.04;
[xt yt n] = hfill2(x,y,nb,lx,hx,nb,lx,hx);
imagesc(xt,yt,n);
set(gca,'YDir','normal');
axis square
xlim([lx hx]);
ylim(xlim());
hold on;
plot([-1 1],[-1 1],'k');
text(lx+(hx-lx)*.05,hx-(hx-lx)*.05,'\deltay (deg)',...
     'HorizontalAlignment','Left','VerticalAlignment','top',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.15,...
     ['r=' sprintf('%0.2f',P(1,2))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.05,...
     ['m=' sprintf('%1.2f',PP(1))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
% Label accordingly
for ii = 1:2
  switch lab{ii}
    case 'bmmeas'
      labstr{ii} = 'Measured Beam Params';
    case 'bmdp'
      labstr{ii} = 'Beam Map Sim dp coeffs';
    case 'realdp'
      labstr{ii} = 'Real dp coeffs';
  end
end
xlabel(labstr{1});
ylabel(labstr{2});

figure(2)
subplot(3,2,5)
scatter(coeff1(:,3),coeff2(:,3),'.');
axis square
xlim([lx hx]);
ylim(xlim());
hold on;
plot([-1 1],[-1 1],'k');
text(lx+(hx-lx)*.05,hx-(hx-lx)*.05,'\deltay (deg)',...
     'HorizontalAlignment','Left','VerticalAlignment','top',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.15,...
     ['r=' sprintf('%0.2f',P(1,2))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.05,...
     ['m=' sprintf('%1.2f',PP(1))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
% Label accordingly
for ii = 1:2
  switch lab{ii}
    case 'bmmeas'
      labstr{ii} = 'Measured Beam Params';
    case 'bmdp'
      labstr{ii} = 'Beam Map Sim dp coeffs';
    case 'realdp'
      labstr{ii} = 'Real dp coeffs';
  end
end
xlabel(labstr{1});
ylabel(labstr{2});
clear x y P PP

% dsig
figure(1)
subplot(3,2,2)
x = cvec(coeff1(:,4));
y = cvec(coeff2(:,4));
cutnan = find(isnan(x)); x(cutnan) = []; y(cutnan) = [];
cutnan = find(isnan(y)); x(cutnan) = []; y(cutnan) = [];
P = corrcoef(x,y);
PP = polyfit(x,y,1);
lx = -0.01;
hx = 0.01;
[xt yt n] = hfill2(x,y,nb,lx,hx,nb,lx,hx);
imagesc(xt,yt,n);
set(gca,'YDir','normal');
axis square
xlim([lx hx]);
ylim(xlim());
hold on;
plot([-1 1],[-1 1],'k');
text(lx+(hx-lx)*.05,hx-(hx-lx)*.05,'\delta\sigma (deg)',...
     'HorizontalAlignment','Left','VerticalAlignment','top',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.15,...
     ['r=' sprintf('%0.2f',P(1,2))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.05,...
     ['m=' sprintf('%1.2f',PP(1))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);

figure(2)
subplot(3,2,2)
scatter(coeff1(:,4),coeff2(:,4),'.');
axis square
xlim([lx hx]);
ylim(xlim());
hold on;
plot([-1 1],[-1 1],'k');
text(lx+(hx-lx)*.05,hx-(hx-lx)*.05,'\delta\sigma (deg)',...
     'HorizontalAlignment','Left','VerticalAlignment','top',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.15,...
     ['r=' sprintf('%0.2f',P(1,2))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.05,...
     ['m=' sprintf('%1.2f',PP(1))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
clear x y P PP

% dp
figure(1)
subplot(3,2,4)
x = cvec(coeff1(:,5));
y = cvec(coeff2(:,5));
cutnan = find(isnan(x)); x(cutnan) = []; y(cutnan) = [];
cutnan = find(isnan(y)); x(cutnan) = []; y(cutnan) = [];
P = corrcoef(x,y);
PP = polyfit(x,y,1);
lx = -0.05;
hx = 0.05;
[xt yt n] = hfill2(x,y,nb,lx,hx,nb,lx,hx);
imagesc(xt,yt,n);
set(gca,'YDir','normal');
axis square
xlim([lx hx]);
ylim(xlim());
hold on;
plot([-1 1],[-1 1],'k');
text(lx+(hx-lx)*.05,hx-(hx-lx)*.05,'\deltap',...
     'HorizontalAlignment','Left','VerticalAlignment','top',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.15,...
     ['r=' sprintf('%0.2f',P(1,2))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.05,...
     ['m=' sprintf('%1.2f',PP(1))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);

figure(2)
subplot(3,2,4)
scatter(coeff1(:,5),coeff2(:,5),'.')
axis square
xlim([lx hx]);
ylim(xlim());
hold on;
plot([-1 1],[-1 1],'k');
text(lx+(hx-lx)*.05,hx-(hx-lx)*.05,'\deltap',...
     'HorizontalAlignment','Left','VerticalAlignment','top',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.15,...
     ['r=' sprintf('%0.2f',P(1,2))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.05,...
     ['m=' sprintf('%1.2f',PP(1))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
clear x y P PP

% dc
figure(1)
subplot(3,2,6)
x = cvec(coeff1(:,6));
y = cvec(coeff2(:,6));
cutnan = find(isnan(x)); x(cutnan) = []; y(cutnan) = [];
cutnan = find(isnan(y)); x(cutnan) = []; y(cutnan) = [];
P = corrcoef(x,y);
PP = polyfit(x,y,1);
lx = -0.05;
hx = 0.05;
[xt yt n] = hfill2(x,y,nb,lx,hx,nb,lx,hx);
imagesc(xt,yt,n);
set(gca,'YDir','normal');
axis square
xlim([lx hx]);
ylim(xlim());
hold on;
plot([-1 1],[-1 1],'k');
text(lx+(hx-lx)*.05,hx-(hx-lx)*.05,'\deltac',...
     'HorizontalAlignment','Left','VerticalAlignment','top',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.15,...
     ['r=' sprintf('%0.2f',P(1,2))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.05,...
     ['m=' sprintf('%1.2f',PP(1))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);

figure(2)
subplot(3,2,6)
scatter(coeff1(:,6),coeff2(:,6),'.');
axis square
xlim([lx hx]);
ylim(xlim());
hold on;
plot([-1 1],[-1 1],'k');
text(lx+(hx-lx)*.05,hx-(hx-lx)*.05,'\deltac',...
     'HorizontalAlignment','Left','VerticalAlignment','top',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.15,...
     ['r=' sprintf('%0.2f',P(1,2))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
text(hx-(hx-lx)*.05,lx+(hx-lx)*.05,...
     ['m=' sprintf('%1.2f',PP(1))],...
     'HorizontalAlignment','Right','VerticalAlignment','bottom',...
     'Color',[0,0,0]);
clear x y P PP

figure(1)
colormap(flipud(gray(256)))
% Label accordingly
switch expt
  case 'bicep2'
    str = ['BICEP2 deprojection coefficients'];
    savename = 'bicep2_';
  case 'keck'
    str1 = ['Keck ' num2str(year) ' '];
    switch perwhat
      case 'perrx'
        str2 = [rxfreq ' '];
      case 'perfreq'
        str2 = [rxfreq ' GHz '];
      otherwise
        str2 = '_';
    end
    str = [str1 str2 ' deprojection coefficients'];
    savename = ['keck' num2str(year) '_' rxfreq];
end
% 'bmmeas','bmdp','realdp'
suffix = [lab{1} '_'  lab{2} '_hist'];
gtitle(str)
print('-depsc2',[savename suffix]);

figure(2)
colormap(flipud(gray(256)))
% Label accordingly
switch expt
  case 'bicep2'
    str = ['BICEP2 deprojection coefficients'];
    savename = 'bicep2_';
  case 'keck'
    str1 = ['Keck ' num2str(year) ' '];
    switch perwhat
      case 'perrx'
        str2 = [rxfreq ' '];
      case 'perfreq'
        str2 = [rxfreq ' GHz '];
      otherwise
        str2 = '_';
    end
    str = [str1 str2 ' deprojection coefficients'];
    savename = ['keck' num2str(year) '_' rxfreq];
end
% 'bmmeas','bmdp','realdp'
suffix = [lab{1} '_'  lab{2} '_scatter'];
gtitle(str)
print('-depsc2',[savename suffix]);


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following is based on Chin Lin's 
% /n/home08/clwong/20140418_beampaper/compdepjcoef_20141229.m
% which probably has way more than necessary!

%{ 
% Make diff beamparams from beamwid file
sig = sqrt((p.fwhm_maj(ind.a).^2 + p.fwhm_min(ind.a).^2 + ...
            p.fwhm_maj(ind.b).^2 + p.fwhm_min(ind.b).^2)/4)/ ...
      (2*sqrt(2*log(2)));
p1.x = p1.r.*cosd(p1.theta);
p1.y = p1.r.*sind(p1.theta);
p1.ukpv = p.ukpv; % except gain, from get_array_info
gain = 1./p1.ukpv;
g_ukpv = (gain(ind.a) - gain(ind.b))./(gain(ind.a) + gain(ind.b))*2;
gain_n(ind.a) = gain(ind.a)./(gain(ind.a) + gain(ind.b))*2;
gain_n(ind.b) = gain(ind.b)./(gain(ind.a) + gain(ind.b))*2;

% Check how different this is from "better" diff params in beamwid file
bmparam(:,1) = g_ukpv;
bmparam(:,2) = (p1.x(ind.a) - p1.x(ind.b));
bmparam(:,3) = (p1.y(ind.a) - p1.y(ind.b));
bmparam(:,4) = p1.sigma(ind.a) - p1.sigma(ind.b);
bmparam(:,5) = p1.p(ind.a) - p1.p(ind.b);
bmparam(:,6) = p1.c(ind.a)-p1.c(ind.b);

% Load real map
load(realmap);

% Average depj coefficients
b = coaddopt.b;
bw = coaddopt.bw;
bw = reshape(bw,[size(b,1),2,1,size(b,4)]);
bw = repmat(bw,[1,1,6,1]);
bnew = squeeze(reshape(b,size(b,1),1,6,[]));
bwnew = squeeze(reshape(bw,size(bw,1),1,6,[]));
for ii = 1:6
  for jj = 1:size(bnew,1)
    V(jj,ii) = nanvar(squeeze(bnew(jj,ii,:)),squeeze(bwnew(jj,ii,:)),1);
  end
end
b = squeeze(nansum(b.*bw,2)./nansum(bw,2));
bw = squeeze(nansum(bw,2));
% Average over time
b = nansum(b.*bw,3)./nansum(bw,3);
depjcoefreal = b;

% Calibrate depj coefs
ukpv = (p.ukpv(ind.a) + p.ukpv(ind.b))/2;
% What the hell is going on here?
switch expt
  case 'bicep2'
    ukpv1=ukpv;
    %ukpv1=3150;
  case 'keck'
    switch year
      case 2012
        ukpv1 = ukpv;
        %ukpv1=3400;
      case 2013
        ukpv1 = ukpv./2900*3400;;
        %ukpv1=3400;
      case 2014
      case 2015
    end
end

% _rc = relgain-calibrated (?)
kk = 1;
depjcoefreal_rc(:,kk) = depjcoefreal(:,kk)*2.*ukpv;
for kk = 2:3;
  depjcoefreal_rc(:,kk) = depjcoefreal(:,kk)*2.*ukpv;
end
kk = 4;
depjcoefreal_rc(:,kk) = depjcoefreal(:,kk)*2.*ukpv./sig;
kk = 5;
depjcoefreal_rc(:,kk) = depjcoefreal(:,kk)*4.*ukpv./(sig.^2);
kk = 6;
depjcoefreal_rc(:,kk) = depjcoefreal(:,kk)*4.*ukpv./(sig.^2);
% Is there another factor of 2 in one of the diff ellips?

% Rotate dx/dy
[THETA,RHO] = cart2pol(depjcoefreal_rc(:,2),depjcoefreal_rc(:,3));
[depjcoefreal_rc(:,2) depjcoefreal_rc(:,3)] = ...
    pol2cart(THETA+p.drumangle(ind.a)*pi/180,RHO);

% Rotate dp/dc 
R(1,1,:) = cosd(2*p.drumangle(ind.a));
R(1,2,:) = -sind(2*p.drumangle(ind.a));
R(2,1,:) = sind(2*p.drumangle(ind.a));
R(2,2,:) = cosd(2*p.drumangle(ind.a));
% Seriously, Chin Lin?
thingold = [depjcoefreal_rc(:,5)'; depjcoefreal_rc(:,6)'];
thingnew = NaN(size(thingold));
for ii = 1:size(thingold,2);
  thingnew(:,ii) = R(:,:,ii)*thingold(:,ii);
end
depjcoefreal_rc(:,5) = thingnew(1,:);
depjcoefreal_rc(:,6) = thingnew(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting time

switch expt
  case 'bicep2'
    [tmp,index{1}] = intersect(ind.a,ind.rgla);         
    expname = {'BICEP2'}                                    
    year = '';
  case 'keck'
    [tmp,index{1}] = intersect(ind.a,ind.rgla);
    for ii = unique(p.rx)'
      [tmp,index{ii+2}] = ...
          intersect(ind.a,ind.rgla(find(p.rx(ind.rgla)==ii)));
    end
    year = num2str(year);
    expname = {'Keck','rx0','rx1','rx2','rx3','rx4'};
    fn = ['Keck ' year]
    texttitle = { fn, [fn ', Rx0'],...
                  [fn ', Rx1'],[fn ', Rx2'],...
                  [fn ', Rx3'],[fn ', Rx4']}
end

Txt = {'relgain depj coef * 2 * ukpv',...
       'dx depj coef * 2 *ukpv / sig',...
       'dy depj coef * 2 *ukpv / sig',...
       'dsig depj coef * 2 *ukpv / sig',...
       'd (+) depj coef * 4 *ukpv / sig^2',...
       'd (x) depj coef * 4 *ukpv / sig^2'};
Txt2 = {'relgain',...
        'dx',...
        'dy',...
        'dsig',...
        'd (+)',...
        'd (x)'};

xtxt = {'\deltag', ... 
        '\deltax (deg)', ...
        '\deltay (deg)',...
        '\delta\sigma (deg)',...
        '\deltap',...
        '\deltac'}

% Plot limits
switch expt
  case 'bicep2'
    A = [-0.05 0.05 -0.05 0.05;...
         -0.02 0.04 -0.02 0.04;...
         -0.02 0.04 -0.02 0.04;...
         -0.01 0.01 -0.01 0.01;...
         -0.05 0.05 -0.05 0.05;...
         -0.05 0.05 -0.05 0.05];
  case 'keck'
    A = [-0.05 0.05 -0.05 0.05;...
         -0.04 0.04 -0.04 0.04;...
         -0.04 0.04 -0.04 0.04;...
         -0.01 0.01 -0.01 0.01;...
         -0.05 0.05 -0.05 0.05;...
         -0.05 0.05 -0.05 0.05];
end

f = figure('visible','off')
clf
setwinsize(gcf,400,600)
switch expt
  case 'bicep2'
    nb = [20 20 20 20 20 20];
  case 'keck'
    nb = [20 25 25 20 20 20];
end
plotparam = [1 4 2 5 3 6];

% Go through each rgl pair (index above)
for jj = 1:length(index)
  clf
  for kk = 1:6
    subplot(3,2,kk)
    mm = plotparam(kk);
    x = bmparam(index{jj},mm);
    y = depjcoefreal_rc(index{jj},mm);
    % Play around before deciding if to cut this
    if makepng
      cutnan = find(isnan(x));
      x(cutnan) = []; y(cutnan) = [];
      cutnan = find(isnan(y));
      x(cutnan) = []; y(cutnan) = [];
      P = corrcoef(x,y);
      PP = polyfit(x,y,1);
    end
    lx = A(mm,1); 
    hx = A(mm,2); 
    ly = A(mm,3); 
    hy = A(mm,4);
    % Actual thing that's plotted
    [xt,yt,n] = hfill2(x,y,nb(mm),lx,hx,nb,ly,hy);
    imagesc(xt,yt,n);
    % Lots of setup just for that...
    set(gca,'YDir','normal')
    hold on
    %      if mm==5
    %	    plot([-1 1],[-1 1]+0.02,'k')
    %      else
    %	    plot([-1 1],[-1 1],'k')
    %      end
    axis square
    xlim([lx,hx]);
    ylim(xlim());
    hold on;
    plot([-1,1],[-1,1],'k');
    hold off;
    if mm < 2 | mm > 3
      tmp=nanmean(dd{mm});
      hold on;
      plot([-1,1],[-1,1] + tmp,'--k');
      hold off;
    end
    text(lx+(hx-lx)*.05,hx-(hx-lx)*.05,xtxt{mm},...
         'HorizontalAlignment','Left','VerticalAlignment','top',...
         'Color',[0,0,0]);
    if makepng
      title({['r=' num2str(P(1,2))],['slope=' num2str(PP(1))]})
    end
    colormap(flipud(gray(256)))
    axis square
    %colorbar
    if kk == 5
      ylabel('deproj. coeff.');
      xlabel('measured');
    end
  end
  
  switch expt
    case 'bicep2'
      annotation('textbox', [0 0.9 1 0.1], ...
                 'String', [expname{jj}], ...
                 'EdgeColor', 'none', ...
                 'HorizontalAlignment', 'center')
      filename = ['realvsbm_bicep2'];
    case 'keck'
      annotation('textbox', [0 0.9 1 0.1], ...
                 'String', [texttitle{jj}], ...
                 'EdgeColor', 'none', ...
                 'HorizontalAlignment', 'center')
      filename = ['realvsbm_' ...
                  experiment year '_' expname{jj}];
  end
  
  if makepng
    mkpng([subdir filename],1)
  else
    print('-depsc2',[subdir filename])
  end

end

%}