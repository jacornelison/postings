% plot_skips('20130509C01_dk293')
% or plot_skips('data/real/20130509C01_dk293_tod.mat')
% or plot_skips(d)
%
% Tool for the reduc czar to examine skips, drops, merge failures, etc.
function plot_skips(tag)

if ischar(tag)
  if exist(tag,'file')
    fname=tag;
    [fdir basename fext]=fileparts(fname);
    tag=strrep(basename,'_tod','');
  else
    fname=fullfile('data','real',tag(1:6),[tag '_tod.mat']);
  end
  load(fname);
else
  liftvars(tag);
  tag='Unknown tag';
end

d.thrs = (d.t-d.t(1))*24;

% Loop over el nods, field scan block
for i=1:3
  switch(i)
    case 1, b.sf=en.sf(1);  b.ef=en.ef(1);  part='Leading elnod';
    case 2, b.sf=fsb.sf(1); b.ef=fsb.ef(1); part='Field scans';
    case 3, b.sf=en.sf(2);  b.ef=en.ef(2);  part='Trailing elnod';
  end

  S = identify_skips(d, b);
  [n,j]=print_skips(S, d, b, tag, part);
  % plot_finder(d,b.sf,b.ef,j,n,tag,part);
  % plot_details(d,b.sf,b.ef,j,n,tag,part);
end

return

%%%%
% Print text output about the identified skips
function [skipnum,skipidx]=print_skips(S, d, fs, tag, part)
  % To report to plotting functions (if enabled)
  skipnum = zeros(1,0);
  skipidx = zeros(1,0);

  [j0,idx] = sort(S.sf);
  j1 = S.ef(idx);
  if length(j0) == 0
    return
  end

  fprintf(1, 'In tag %s, %s:\n', tag, part);

  hrsblock = (d.t(fs.sf)-d.t(1))*24 + ...
      (d.t(fs.ef)-d.t(fs.sf))*24*linspace(0,1,fs.ef-fs.sf+1);

  % To keep track of unique events
  lastskip = 0;
  onskip = 0;
  for ii=1:length(j0)
    % Number common skips together by keeping track of the time at which each
    % problem is identified. Update as a later time is reached.
    if j0(ii) > lastskip
      onskip = onskip + 1;
      lastskip = j0(ii);
      % Plotting uses block-relative indexing, so fix up when appending to
      % skipidx.
      skipnum = [skipnum, onskip];
      skipidx = [skipidx, lastskip-fs.sf+1];
    end

    % Attempt to find a valid time index, backing up from the identified fault
    % until a non-zero time sample is identified.
    jj = j0(ii);
    while d.t(jj) == 0 && jj > 0
      jj = jj - 1;
    end
    % If we've baked up beyond the beginning of the block, then take the
    % first valid time sample we can find.
    if jj < fs.sf
      jj = find(d.t(fs.sf:fs.ef)~=0, 1);
    end

    % Get a printable time string.
    [YY MM DD HH NN SS] = mjd2date(d.t(jj));
    dn = datenum(YY, MM, DD, HH, NN, SS);
    dstr = datestr(dn, 'yyyy-mmm-dd:HH:MM:SS');
    % Then print a message to the screen
    fprintf(1, '    (%d) %s: %s (t=%0.4f hrs; %s)\n', onskip, ...
        S.name{idx(ii)}, S.desc{idx(ii)}, hrsblock(j0(ii)-fs.sf+1), dstr);
  end
  disp(' ');
return

%%%%
% Make a zoomed-out finder plot
function plot_finder(d,sf,ef,j,n,tag,part)
feat=nrep(d.array.frame.features,20);
feat=nanmedian(feat(sf:ef));
iselnod=(feat==9);
if iselnod
  pos=d.pointing.hor.el(sf:ef);
else
  pos=d.pointing.hor.az(sf:ef);
end
tsmooth=d.thrs(sf)+(d.thrs(ef)-d.thrs(sf))/(ef-sf)*(0:(ef-sf));
% tsmooth=sf:ef;
figure;
plot(tsmooth,pos,'r-');
grid on;
hold on;
axax=axis;
ylim([axax(3) axax(4)+(axax(4)-axax(3))*0.25]);
axax=axis;
for i=1:length(j)
  ypos=axax(3)+(axax(4)-axax(3))*(0.9-0.05*mod(n(i)-1,4));
  text(tsmooth(j(i)),ypos,[num2str(n(i))],'horizontalalignment','center','verticalalignment','bottom');
end
title([strrep(tag,'_','\_') ', ' part]);
% xlabel('sample number');
xlabel('time / hr');
return

%%%%
% Make a zoomed-in plot for each skip
function plot_details(d,sf,ef,j,n,tag,part)
feat=nrep(d.array.frame.features,20);
feat=nanmedian(feat(sf:ef));
iselnod=(feat==9);
if iselnod
  pos=d.pointing.hor.el(sf:ef);
  labl='el / deg';
else
  pos=d.pointing.hor.az(sf:ef);
  labl='az / deg';
end
tsmooth=d.thrs(sf)+(d.thrs(ef)-d.thrs(sf))/(ef-sf)*(0:(ef-sf));
nmce=size(d.mce0.syncBox.sampleNumber,2);
jdx=1:(ef-sf);
jdx2=1:length(d.t);
for i=1:length(j)
  figure;
  cc=abs(jdx-j(i))<100;
  cc2=abs(jdx2-sf-j(i))<100;
  strip_plot_helper(tsmooth(cc),pos(cc),1,nmce+3,labl);
  title([strrep(tag,'_','\_') ', ' part ' skip #' num2str(i)]);
  strip_plot_helper(tsmooth(cc),d.antenna0.time.utcfast(cc2),2,nmce+3,'UTC');
  strip_plot_helper(tsmooth(cc),d.antenna0.syncBox.sampleNumber(cc2),3,nmce+3,'ant sync');
  for rx=1:nmce
    strip_plot_helper(tsmooth(cc),d.mce0.syncBox.sampleNumber(cc2,rx),3+rx,nmce+3,['mce' num2str(rx-1) ' sync']);
  end
end
return

%%%%%
function strip_plot_helper(x,y,i,n,labl)
axstrip(n,i);
plot(x,y,'r-');
xlim([min(x) max(x)]);
yl=[min(y(y~=0)) max(y(y~=0))];
if any(y==0)
  if yl(1)>0
    yl(1)=2*yl(1)-yl(2);
  else
    yl(2)=2*yl(2)-yl(1);
  end
end
ylim(yl);
text(0.1,0.9,labl,'units','normalized');
set(gca,'yticklabel','');
if i<n
  set(gca,'xticklabel','');
else
  xlabel('time / hr');
end
return

%%%%%
% Remaining functions below are nuts & bolts utilities
%

%%%%
function ax=axstrip(m,i)
% Fractional border at LRBT
b=[0.1 0.06 0.1 0.07];
ax=axes('position',[b(1) b(3)+(1-b(3)-b(4))/m*(m-i) 1-b(1)-b(2) (1-b(3)-b(4))/m]);
return
