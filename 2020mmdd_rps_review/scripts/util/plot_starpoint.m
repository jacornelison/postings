function plot_starpoint(csvfile,outdir,onlyhtml)
% plot_starpoint(csvfile,outdir,onlyhtml)
%
% Make starpointing plots and appropriate html for all pointing model data.
% 
% ex. plot_starpoint('aux_data/pointingmodels/pointingmodel_complete_v1.csv','reducplots/starpoint',0);

if ~exist('csvfile','var')
  csvfile=[];
end
if ~exist('outdir','var')
  outdir=[];
end
if ~exist('onlyhtml','var')
  onlyhtml=0;
end

if isempty(csvfile)
  csvfile='aux_data/pointingmodels/pointingmodel_complete_v1.csv';
end
if isempty(outdir)
  outdir='starpoint_plots';
end

[p,k]=ParameterRead(csvfile);
dn=datenum(p.meas_name,'yymmdd HH:MM:SS');

if ~onlyhtml
  figure(1)
  clf
  setwinsize(gcf,1000,1100);
  
  % axis tilts
  subplot(4,1,1);
  hold on;
  plot(dn,p.az_tilt_ha,'.-');plot(dn,p.az_tilt_lat','.-g');plot(dn,p.el_tilt,'.-r');
  hold off;
  datetick('x','mmm/yy');grid on;
  ax0=gca;
  title('axis tilts (blue:az\_tilt\_ha, green:az\_tilt\_lat, red:el\_tilt');
  ylabel('tilt');
  
  % collim mag/dir
  subplot(4,1,2);
  % make collim dir >0<180
  ind=find(p.collim_dir>180);
  p.collim_dir(ind)=p.collim_dir(ind)-180;
  p.collim_mag(ind)=-p.collim_mag(ind);
  ax1=plottwo(dn,p.collim_mag,p.collim_dir,'collim_mag','collim_dir',ax0);
  title('collimation');
  
  % encoder zeros
  subplot(4,1,3)
  ax2=plottwo(dn,p.az_zero,p.el_zero,'az_zero','el_zero',ax0);
  title('encoder zero points');
  
  % Nstars and residuals
  subplot(4,1,4)
  ax3=plottwo(dn,p.rms_res,p.numpts,'rms_res','Nstars',ax0);
  title('rms residuals / Nstars');
  
  % now loop over dates, adding a dashed line for each dates
  for i=1:numel(dn)
    
    clear h
    
    set(1,'CurrentAxes',ax0);
    yl=get(gca,'Ylim');
    hold on;h(1)=plot(ax0,[dn(i),dn(i)],[yl(1),yl(2)],'--k');hold off;
    
    set(1,'CurrentAxes',ax1(1));
    yl=get(gca,'Ylim');
    hold on;h(2)=plot(gca,[dn(i),dn(i)],[yl(1),yl(2)],'--k');hold off;
    
    set(1,'CurrentAxes',ax2(1));
    yl=get(gca,'Ylim');
    hold on;h(3)=plot(gca,[dn(i),dn(i)],[yl(1),yl(2)],'--k');hold off;
    
    set(1,'CurrentAxes',ax3(1));
    yl=get(gca,'Ylim');
    hold on;h(4)=plot(gca,[dn(i),dn(i)],[yl(1),yl(2)],'--k');hold off;
    
    plotdate=datestr(dn(i),'yymmddHHMMSS');
    plottit=sprintf('%s/%s_all.png',outdir,plotdate);
    mkpng(plottit);
    setpermissions(plottit);
    
    for k=1:numel(h)
      delete(h(k));
    end
    
  end

end % ~only html


% make browser html
make_starpoint_html(dn,outdir);


return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ax=plottwo(x,y1,y2,ylab1,ylab2,ax0)

[ax,h1,h2]=plotyy(x,y1,x,y2);
set(h1,'Marker','.')
set(h2,'Marker','.')
datetick(ax(1),'x','mmm/yy');
datetick(ax(2),'x','mmm/yy');
xl=get(ax0,'XLim');
xlim(ax(1),xl);
xlim(ax(2),xl);
set(ax(1),'XTick',[]);
set(ax(2),'XTick',[]);
set(ax(1),'XTick',get(ax0,'XTick'));
set(ax(1),'XTickLabel',get(ax0,'XTickLabel'));
grid on;
set(get(ax(1),'Ylabel'),'String',ylab1)
set(get(ax(2),'Ylabel'),'String',ylab2)
set(get(ax(1),'Ylabel'),'Interpreter','none')
set(get(ax(2),'Ylabel'),'Interpreter','none')

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_starpoint_html(dn,outdir)

% make index
fname=sprintf('%s/index.html',outdir);
h=fopen(fname,'w');

fprintf(h,'<SCRIPT LANGUAGE="JavaScript">\n');
fprintf(h,'<!--\n\n');

dt=datestr(dn(end),'yymmddHHMMSS');
fprintf(h,'date=''%s'';\n',dt);
fprintf(h,'plottype=''all'';\n\n');

fprintf(h,'fig_fname=date+''_''+plottype+''.png'';\n\n');

fprintf(h,'function plupdate(){\n');
fprintf(h,'  fig_fname=date+''_''+plottype+''.png'';\n');
fprintf(h,'  plotpage.document["fig"].src=fig_fname;\n');
fprintf(h,'}\n');
fprintf(h,'function set_date(xx){\n');
fprintf(h,'  date=xx;\n');
fprintf(h,'  plupdate();\n');
fprintf(h,'}\n');
fprintf(h,'function set_plottype(xx){\n');
fprintf(h,'  plottype=xx;\n');
fprintf(h,'  plupdate();\n');
fprintf(h,'}\n\n');

fprintf(h,'//-->\n');
fprintf(h,'</SCRIPT>\n\n');
 
fprintf(h,'<html>\n\n');

fprintf(h,'<head><title>Starpointing</title></head>\n\n');

fprintf(h,'<frameset noresize="noresize" cols="200,*">\n\n');

fprintf(h,'<frame src="starpointing_dates.html" name="dates">\n');
fprintf(h,'<frame src="starpointing_plots.html" name="plotpage">\n\n');

fprintf(h,'</frameset>\n\n');

fprintf(h,'</html>\n');

fclose(h);



% make dates panel
fname=sprintf('%s/starpointing_dates.html',outdir);
h=fopen(fname,'w');

fprintf(h,'<SCRIPT LANGUAGE="JavaScript">\n');
fprintf(h,'<!--\n\n');

fprintf(h,'function set_date(date){\n');
fprintf(h,'  parent.set_date(date);\n');
fprintf(h,'}\n');
fprintf(h,'//-->\n');
fprintf(h,'</SCRIPT>\n\n');

fprintf(h,'<html>\n');
fprintf(h,'<body bgcolor="#d0d0d0">\n');
fprintf(h,'<pre>\n');
fprintf(h,'Starpointing dates:<hr>\n');
for i=1:numel(dn)
  s1=datestr(dn(end-i+1),'yymmddHHMMSS');
  s2=datestr(dn(end-i+1),'yyyymmdd HHMMSS');
  fprintf(h,'<a href="javascript:set_date(''%s'');">%s</a>\n',s1,s2);
end

fprintf(h,'</pre>\n\n');

fprintf(h,'</html>\n');
fclose(h);



% make plots page
fname=sprintf('%s/starpointing_plots.html',outdir);
h=fopen(fname,'w');

fprintf(h,'<SCRIPT LANGUAGE="JavaScript">\n');
fprintf(h,'<!--\n\n');

fprintf(h,'function set_plottype(plottype){\n');
fprintf(h,'  parent.set_plottype(plottype)\n');
fprintf(h,'}\n');
fprintf(h,'//-->\n');
fprintf(h,'</SCRIPT>\n\n');

fprintf(h,'<html>\n\n');
	
fprintf(h,'<h2><center><b>Starpointing</b></center></h2></td>\n\n');

fprintf(h,'<center>\n');
fprintf(h,'<a href="javascript:set_plottype(''all'');">All</a> |\n');
fprintf(h,'<a href="javascript:set_plottype(''1'');">online</a> |\n');
fprintf(h,'<a href="javascript:set_plottype(''2'');">re-fit</a> |\n');
fprintf(h,'<a href="javascript:set_plottype(''3'');">exclude outliers</a>\n');
fprintf(h,'</center>\n\n');

fprintf(h,'<p>\n\n');

dt=datestr(dn(end),'yymmddHHMMSS');
fprintf(h,'<img src="%s_all.png" name="fig">\n\n',dt);

fprintf(h,'<SCRIPT LANGUAGE="JavaScript">\n');
fprintf(h,'<!--\n');
fprintf(h,'parent.plupdate();\n');
fprintf(h,'//-->\n');
fprintf(h,'</SCRIPT>\n\n');

fprintf(h,'</html>\n');

fclose(h);

return
