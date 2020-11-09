function make_cryo_html(dayend)
plotdir='cryoplots/';
htmldir='cryoplots/';

tstart=datenum('20120212','yyyymmdd');
tend=datenum(num2str(dayend),'yyyymmdd');

tlist=tstart:7:tend;
tlist=datestr(tlist,'yyyymmdd');
[mm nn]=size(tlist);

% List each cycle and in-between
page_list=[];
for ii=1:mm
  if exist(['cryoplots/' tlist(ii,:) '.png'])
    page_list(end+1).name=tlist(ii,:);
    page_list(end).str=tlist(ii,:);
  end
end

% Each individual html page
for ii=1:length(page_list)
  h=fopen(fullfile(htmldir,[page_list(ii).name '.html']),'wt');
  fprintf(h,'<html>\n\n<body>\n\n');
  fprintf(h,'    <center><h3>Cryo data</h3></center>\n');
  fprintf(h,'    <center><h3>');
  if ii>1
    fprintf(h,'[<a href="%s.html">prev</a>]',page_list(ii-1).name);
  else
    fprintf(h,'<b>[prev]</b>');
  end
  fprintf(h,' [%s] ',page_list(ii).str);
  if ii<length(page_list)
    fprintf(h,'[<a href="%s.html">next</a>]',page_list(ii+1).name);
  else
    fprintf(h,'<b>[next]</b>');
  end
  fprintf(h,'</center>\n');
  fprintf(h,'<center><h3>');
  fprintf(h,'</center>\n');
  fprintf(h,'<table border="0" cellspacing="0" cellpadding="0">\n');
  fprintf(h,'<tr>\n');
  plotlist=dir(fullfile(plotdir,[page_list(ii).name '.png']));
  for jj=1:length(plotlist)
    fprintf(h,'  <td>\n');
    fprintf(h,'    <a href="%s" target="_blank">',plotlist(jj).name);
    fprintf(h,'<img src="%s">\n',plotlist(jj).name);
    fprintf(h,'  </td>\n');
  end
  fprintf(h,'</tr>\n');
  fprintf(h,'</table>\n\n</body>\n\n</html>\n\n');
  fclose(h);
end

% Left pane with list of cycles and examiner pages
h=fopen(fullfile(htmldir,'plot_index.html'),'wt');
fprintf(h,'<html>\n\n');
fprintf(h,'<body bgcolor="#d0d0d0">\n\n');
fprintf(h,'Go to:\n');
fprintf(h,'<pre>\n');
fprintf(h,'<b><font color="blue">Cryo data</font></b>\n');
fprintf(h,'<a href="20120213.html" target="plotpage"><font color="blue">Cumulative</a></a>\n');
for ii=length(page_list):-1:1
  fprintf(h,'<a href="%s.html" target="plotpage"><font color="blue">%s</a></a>\n',page_list(ii).name,page_list(ii).str);
end
fprintf(h,'</pre>\n\n');
fprintf(h,'</body>\n\n');
fprintf(h,'</html>\n\n');
fclose(h);

% Top-level index page
h=fopen(fullfile(htmldir,'index.html'),'wt');
fprintf(h,'<html>\n\n');
fprintf(h,'<frameset noresize="noresize" cols="200,*">\n\n');
fprintf(h,'<frame src="plot_index.html">\n');
fprintf(h,'<frame src=20120213.html name="plotpage">\n\n');
fprintf(h,'</frameset>\n\n</html>\n\n');
fclose(h);

return
