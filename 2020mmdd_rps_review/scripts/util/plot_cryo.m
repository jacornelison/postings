function plot_cryo(daystart,dayend)
ofint=daystart:1:dayend;
ndays=length(datenum(num2str(daystart),'yyyymmdd'):datenum(num2str(dayend),'yyyymmdd'));

% Spacing of x-ticks
if length(ofint)<10
  sp=0.5;
else
  sp=ndays/14;
end

% Get cold head data:
cryo=dir('cryo/201*cryo.csv');
cryonames=[{cryo(:).name}];

% Get HDI data:
hdi=dir('cryo/201*hdi.csv');
hdinames=[{hdi(:).name}];

% Get pressure data:
p=dir('cryo/201*press.csv');
pnames=[{p(:).name}];

% Get flow meter data:
fl=dir('cryo/201*flow.csv');
fnames=[{fl(:).name}];

% Just grab the data in the date range of interest
for ii=1:length(cryo)
  dates(ii)=str2num(cryonames{ii}(1:8));
end
cryo=cryo(ismember(dates,ofint));
cryo=cryo([cryo.bytes]>166);

dates=[];
for ii=1:length(hdi)
  dates(ii)=str2num(hdinames{ii}(1:8));
end
hdi=hdi(ismember(dates,ofint));
dates_hdi=dates(ismember(dates,ofint));

dates=[];
for ii=1:length(p)
  dates(ii)=str2num(pnames{ii}(1:8));
end
p=p(ismember(dates,ofint));

dates=[];
for ii=1:length(fl)
  dates(ii)=str2num(fnames{ii}(1:8));
end
fl=fl(ismember(dates,ofint));

% Parse the cryo data
[t tnum c titles]=parse_cryo(cryo);

% Create figure
figure1 = figure;
clf

% Create subplot
subplot1 = subplot_stack(3,1,1,figure1);
% Uncomment the following line to preserve the X-limits of the axes
box(subplot1,'on');
hold(subplot1,'all');


% Create multiple lines using matrix input to plot
plot1 = plot(tnum,c,'Parent',subplot1,'Marker','.','LineStyle','none');

for ii=1:length(plot1)
  set(plot1(ii),'DisplayName',titles{ii});
end

tstart=datenum(num2str(daystart),'yyyymmdd');
tend=datenum(num2str(dayend),'yyyymmdd');
xt=tstart:sp:tend;
grid on
t=datestr(xt,6);
set(gca,'XTick',xt); set(gca,'XTickLabel',t);

title('Cold head temp [K], HDI Level [L], and Tank Gauge Pressure [psi]','interpreter','none')
xlim([tstart (tend)]) 
ylim([3 8])
ylabel('T [K]')

hold(subplot1,'off')



if ~isempty(p)
  % Parse pressure data
  [t tnum c titles]=parse_p(p);
  [t2 flow titles2]=parse_flow(fl);
  flow(abs(flow)>100) = NaN;

  % Create subplot
  subplot3 = subplot_stack(3,1,3,figure1);
  
  if ~isempty(titles)
    box(subplot3,'on');
    hold(subplot3,'all');

    atm=c(strcmpi(titles,'Atmosphere/psi'),:);
    atmav = nanmean(atm);
    atm=repmat(atm,3,1);

    c=c(~strcmpi(titles,'Atmosphere/psi'),:);
    titles=titles(~strcmpi(titles,'Atmosphere/psi'));
    c_rel=c-atm;

    % Reorder to match other plots:
    c2(1,:)=c_rel(strcmpi(titles,'Theo/psi'),:);
    c2(2,:)=c_rel(strcmpi(titles,'Alvin/psi'),:);
    c2(3,:)=c_rel(strcmpi(titles,'Simon/psi'),:);
    titles={'Theodore/psi','Alvin/psi','Simon/psi'};

    % Create plot
    plot3=plot(tnum,c2,'Parent',subplot3,'Marker','.','LineStyle','none');
    plot4=plot(t2,flow(1,:)-atmav,'b');
    plot5=plot(t2,flow(2,:)-atmav,'Color',[0 .5 0]);
    plot6=plot(t2,flow(3,:)-atmav,'r');
    set(plot3(1),'DisplayName',titles{1});
    set(plot3(2),'DisplayName',titles{2});
    set(plot3(3),'DisplayName',titles{3});
    set(plot4,'DisplayName',[titles{1} ' (flow meter)']);
    set(plot5,'DisplayName',[titles{2} ' (flow meter)']);
    set(plot6,'DisplayName',[titles{3} ' (flow meter)']);
    
    xlim([tstart (tend)]) 
    ylim([0 13])
    t=datestr(xt,6);num2str(daystart)
    set(gca,'XTick',xt); set(gca,'XTickLabel',t);
    grid on

    % Create ylabel
    ylabel('Gauge Pressure [psi]');

    % Create legends
    legend(subplot3,'show');
  end
end
legend(subplot1,'show');
hold(subplot3,'off')

% Create subplot
subplot2 = subplot_stack(3,1,2,figure1);
box(subplot2,'on');
hold(subplot2,'all');

% If no hdi data exists, just jump to the plot generation
if ~isempty(hdi)
  change_date = 20120213;
  if ~isempty(dates_hdi(dates_hdi<change_date));
    [tnum2 h titles]=parse_hdi(hdi(dates_hdi<change_date));
  else
    tnum2=[]; h=[];
  end
  if ~isempty(dates_hdi(dates_hdi>=change_date));
    [dum tnum3 h2 titles]=parse_cryo(hdi(dates_hdi>=change_date));
  
    % Reorder to match other plots:
    h3(1,:)=h2(strcmpi(titles,'I_Theodore/mm'),:);
    h3(2,:)=h2(strcmpi(titles,'I_Alvin/mm'),:);
    h3(3,:)=h2(strcmpi(titles,'I_Simon/mm'),:);
  else
    tnum3=[]; h3=[];
  end
  tnum2=[tnum2 tnum3];
  h = [h h3];
  titles={'Theodore/L','Alvin/L','Simon/L'};

  for ii=1:3
    v(ii,:)=calc_tank_vol(h(ii,:));
    cc = t2~=0;
    flowint(ii,:)=interp1(t2(cc),flow(ii,cc),tnum2);
  end

  % for the data for which the flow meter is missing, just put in 8psi relative, 17.5 absolute:
  flowint(isnan(flowint)) = 17.5; 
  
  % Calculate pressure correction:
  c4 = 1.51616E-11;
  c3 = -8.705869E-9;
  c2 = 1.483225E-6;
  c1 = -3.004364E-4;
  c0 = 0.1476713;
  m1 = 7.595161;

  vu = v; %uncorrected volume
  flowint = flowint*6.894; % psi to kPa
  v = v.*m1.*(c0+c1*flowint+c2*flowint.^2+c3*flowint.^3+c4*flowint.^4);
  
  v1 = v(1,~isnan(v(1,:))); v2 = v(2,~isnan(v(2,:))); v3 = v(3,~isnan(v(3,:)));
  if isempty(v1), v1=0; else usage(1) = 0; end
  if isempty(v2), v2=0; else usage(2) = 0; end
  if isempty(v3), v3=0; else usage(3) = 0; end
  
  usage(1) = v1(end)-v1(1);
  usage(2) = v2(end)-v2(1);
  usage(3) = v3(end)-v3(1);
  
  ttemp = tnum2(~isnan(v(1,:)));
  
  rate = usage./(ttemp(end)-ttemp(1));
  totalv = v1(end)+v2(end)+v3(end);
  
  % Calculate usage:
%   usage=(v(:,end)-v(:,1));
%   rate=(vn(:,1)-vn(:,end))./(tnum2(1)-tnum2(end));
%   % And total volume:
%   totalv=sum(vn(:,end));
  
  % Create plot
  plot2=plot(tnum2,v,'Parent',subplot2,'Marker','.','LineStyle','none');
  plot4=plot(tnum2,vu(1,:),'ob');
  plot5=plot(tnum2,vu(2,:),'o','Color',[0 .5 0]);
  plot6=plot(tnum2,vu(3,:),'or');
  set(plot4,'DisplayName',[titles{1} ' (no P correction)']);
  set(plot5,'DisplayName',[titles{2} ' (no P correction)']);
  set(plot6,'DisplayName',[titles{3} ' (no P correction)']);
  set(plot2(1),'DisplayName',titles{1});
  set(plot2(2),'DisplayName',titles{2});
  set(plot2(3),'DisplayName',titles{3});
  
  text(tstart+sp,2700,['Theo change: ' num2str(usage(1),'%+6.2f') ' L'],'fontsize',14)
  text(tstart+sp,2550,['Alvin change: ' num2str(usage(2),'%+6.2f') ' L'],'fontsize',14)
  text(tstart+sp,2400,['Simon change: ' num2str(usage(3),'%+6.2f') ' L'],'fontsize',14)
  text(tstart+sp,2250,['Total change: ' num2str(sum(usage),'%+6.2f') ' L'],'fontsize',14)
  text(tstart+sp,2100,['Total Volume: ' num2str(totalv,'%6.0f') ' L'],'fontsize',14)
  
  text(tstart+4*sp,2700,['Theo rate: ' num2str(rate(1),'%+6.2f') ' L/day'],'fontsize',14)
  text(tstart+4*sp,2550,['Alvin rate: ' num2str(rate(2),'%+6.2f') ' L/day'],'fontsize',14)
  text(tstart+4*sp,2400,['Simon rate: ' num2str(rate(3),'%+6.2f') ' L/day'],'fontsize',14)
  text(tstart+4*sp,2250,['Total rate: ' num2str(sum(rate),'%+6.2f') ' L/day'],'fontsize',14)
  
  xlim([tstart (tend)]) 
  ylim([0 4000])
  t=datestr(xt,6);num2str(daystart)
  set(gca,'XTick',xt); set(gca,'XTickLabel',t);
  grid on

  % Create ylabel
  ylabel('Volume [L]');

  % Create legends
  legend(subplot2,'show');
end
hold(subplot2,'off')

setwinsize(figure1,1600,1000);
mkpng(['cryoplots/' num2str(daystart)],1);


function [t tnum c titles]=parse_cryo(cryo)
c=[];
hdr=[];
t=[];
tnum=[];
for ii=1:length(cryo)
  clear tnumd
  try
    [d vars]=importfile(['cryo/' cryo(ii).name]);
  catch
    continue
  end
 cd=d.data;
 hdrd=d.textdata(8,:);
 td=d.textdata(9:end,1);
 for jj=1:length(td)
   tnumd(jj)=datenum(td{jj},'yyyy-mm-dd HH:MM:SS');
 end
 if length(cd)~=length(tnumd)
   continue
 end
 try
   c=vertcat(c,cd);
 catch
   continue
 end
 hdr=vertcat(hdr,hdrd);
 t=vertcat(t,td);
 tnum=horzcat(tnum,tnumd);
end
c=c';
t=t';

titles=hdr(end,2:4);
titles{1}=titles{1}(3:end);
titles{2}=titles{2}(3:end);
titles{3}=titles{3}(3:end);

function [t h titles]=parse_hdi(hdi)
h=[];
t=[];
tnum=[];
keep=[];
i=1;
for ii=1:length(hdi)
  fid=fopen(['cryo/' hdi(ii).name]);
  if datenum(hdi(ii).name(1:8),'yyyymmdd') >= 734827
    fgetl(fid);  fgetl(fid);  fgetl(fid);
    titles{1}=fgetl(fid); titles{1}=titles{1}(17:21);
    if strcmp(titles{1},'mon/p')
      continue
    end
    titles{2}=fgetl(fid); titles{2}=titles{2}(17:21);
    titles{3}=fgetl(fid); titles{3}=titles{3}(17:21);
    fgetl(fid);  fgetl(fid);
    while 1
      l=fgetl(fid);
      if ~ischar(l), break, end    
      try
      if strmatch(l(1:2),'20')
        t(i)=datenum(l(1:19),'yyyy-mm-dd HH:MM:SS');
        if strcmpi(titles{1},'simon') && strcmpi(titles{2},'alvin') && strcmpi(titles{3},'theod')
          ind=strfind(l,'mm');
          c1(i)=str2double(l(ind-4:ind-1));
          l=fgetl(fid);
          ind=strfind(l,'mm');
          c2(i)=str2double(l(ind-4:ind-1));
          l=fgetl(fid);
          if length(l)>2
            ind=strfind(l,'mm');
            c3(i)=str2double(l(ind-4:ind-1)); 
            fgetl(fid);
          else
            c3(i)=NaN;
          end
        end
      end
      catch
        warning(['There is a problem with ' hdi(ii).name '.  Skipping for now. You may want to inspect file by hand.'])
        c1(i)=NaN; c2(i)=NaN; c3(i)=NaN;
      end
      i=i+1;
    end
  else
    l=fgetl(fid);
    while ischar(l)
      l=fgetl(fid);
      if strmatch('# Column 2',l)
        name=l;
      end
      if strmatch('201',l)
        t(i)=datenum(l(1:19),'yyyy-mm-dd HH:MM:SS');
        if length(l)>23
          if ~isempty(strfind(name,'imon'))
            c1(i)=str2double(l(23:26)); c2(i)=NaN; c3(i)=NaN;
          elseif ~isempty(strfind(name,'lvin'))
            c2(i)=str2double(l(23:26)); c1(i)=NaN; c3(i)=NaN;
          elseif ~isempty(strfind(name,'heodore')) 
            c3(i)=str2double(l(23:26)); c1(i)=NaN; c2(i)=NaN;
          end
        else
          c1(i)=NaN; c2(i)=NaN; c3(i)=NaN;
        end
        i=i+1;
      end
    end
  end
end
i=0;

h=[c3; c2; c1];

titles={'Theodore/L','Alvin/L','Simon/L'};


function [t tnum p titles]=parse_p(pv)
p=[];
hdr=[];
t=[];
tnum=[];
for ii=1:length(pv)
  clear tnumd
  try
    [d vars]=importpress(['cryo/' pv(ii).name]);
  catch
    continue
  end
 pd=d.data;
 hdrd=d.textdata(9,:);
 td=d.textdata(10:end,1);
 for jj=1:length(td)
   tnumd(jj)=datenum(td{jj},'yyyy-mm-dd HH:MM:SS');
 end
 if length(pd)~=length(tnumd)
   continue
 end
 try
   p=vertcat(p,pd);
 catch
   continue
 end
 hdr=vertcat(hdr,hdrd);
 t=vertcat(t,td);
 tnum=horzcat(tnum,tnumd);
end
p=p';
t=t';

if ~isempty(p)
  titles=hdr(end,2:5);
  titles{1}=titles{1}(3:end);
  titles{2}=titles{2}(3:end);
  titles{3}=titles{3}(3:end);
  titles{4}=titles{4}(3:end);
else
  titles=[];
end



function [t flow titles]=parse_flow(pv)
p=[];
t=[];
i=1;

for ii=1:length(pv)
  fid=fopen(['cryo/' pv(ii).name]);
  fgetl(fid);  fgetl(fid);  fgetl(fid);
  titles{1}=fgetl(fid); titles{1}=titles{1}(15:19);
  titles{2}=fgetl(fid); titles{2}=titles{2}(15:19);
  titles{3}=fgetl(fid); titles{3}=titles{3}(15:19);
  
  fgetl(fid);  fgetl(fid);
  while 1
    l=fgetl(fid);
    if ~ischar(l), break, end    
    if strmatch(l(1:2),'20')
      t(i)=datenum(l(1:19),'yyyy-mm-dd HH:MM:SS');
      ind=strfind(l,'T');
      c1(i)=str2double(l(ind+2:ind+9));
      ind=strfind(l,'A');
      c2(i)=str2double(l(ind+2:ind+9));
      ind=strfind(l,'S');
      c3(i)=str2double(l(ind+2:ind+9));
    end
    i = i+1;
  end
end

flow=[c1; c2; c3];

titles={'Theodore/psi','Alvin/psi','Simon/psi'};

if isempty(p)
  titles=[];
end

function [d vars]=importfile(fileToRead1)
%IMPORTFILE(FILETOREAD1)
%  Imports data from the specified file
%  FILETOREAD1:  file to read

%  Auto-generated by MATLAB on 10-Nov-2011 23:52:27

DELIMITER = ',';
HEADERLINES = 8;

% Import the file
d = importdata(fileToRead1, DELIMITER, HEADERLINES);

% Create new variables in the base workspace from those fields.
vars = fieldnames(d);
for i = 1:length(vars)
  assignin('base', vars{i}, d.(vars{i}));
end


function [d vars]=importpress(fileToRead1)
%IMPORTFILE(FILETOREAD1)
%  Imports data from the specified file
%  FILETOREAD1:  file to read

%  Auto-generated by MATLAB on 10-Nov-2011 23:52:27

DELIMITER = ',';
HEADERLINES = 9;

% Import the file
d = importdata(fileToRead1, DELIMITER, HEADERLINES);

% Create new variables in the base workspace from those fields.
vars = fieldnames(d);
for i = 1:length(vars)
  assignin('base', vars{i}, d.(vars{i}));
end


function [d vars]=importflow(fileToRead1)
%IMPORTFILE(FILETOREAD1)
%  Imports data from the specified file
%  FILETOREAD1:  file to read

%  Auto-generated by MATLAB on 10-Nov-2011 23:52:27

DELIMITER = ',';
HEADERLINES = 8;

% Import the file
d = importdata(fileToRead1, DELIMITER, HEADERLINES);

% Create new variables in the base workspace from those fields.
vars = fieldnames(d);
for i = 1:length(vars)
  assignin('base', vars{i}, d.(vars{i}));
end
