% write_tag_html(pt,ps,ks,tag,html_dir,file_url,experiment,do_show_all)
%
% Writes the tag browser HTML page for a single tag.
%
%   pt is the structure from the tag list file.
%   [ps,ks] are the structures from the tag stats file.
%   tag is either a tag name, or the index of a tag within pt.tag.
%   html_dir is the directory in which to write the file.
%   file_url is the URL at which the reducplots live, for image links.
%   experiment is bicep2, keck, or bicep3.  Controls per-RX or per-MCE plots.
%   do_show_all is true if per-phase cals, cryo, etc. are included
%     in "prev phase", "next tag", etc. links.
function write_tag_html(pt,ps,ks,tag,html_dir,file_url,experiment,do_show_all)

if ischar(tag)
  ii=strmatch(tag,pt.tag,'exact');
  if isempty(ii)
    error(['Tag ' tag ' not found in tag list.']);
  end
else
  ii=tag;
end

jj=strmatch(pt.tag{ii},ps.tag,'exact');
if isempty(jj)
  jj=NaN;
else
  jj=jj(1);
end

if ~isfield(pt,'good')
  pt.good{ii}='';
end

on_month=pt.tag{ii}(1:6);

if ~exist(fullfile(html_dir,on_month),'dir')
  mkdir(fullfile(html_dir,on_month));
  fileattrib(fullfile(html_dir,on_month),'+w','g');
end
if ~exist(fullfile(html_dir,on_month,'browser'),'dir')
  mkdir(fullfile(html_dir,on_month,'browser'));
  fileattrib(fullfile(html_dir,on_month,'browser'),'+w','g');
end

htmlfname=fullfile(html_dir,on_month,'browser',[pt.tag{ii} '.html']);

h=fopen(htmlfname,'wt');
if h==-1
  error(['Failed to open file ' htmlfname]);
end

fprintf(h,'<html>\n\n');
fprintf(h,'<head>\n\n');
fprintf(h,'</head>\n<body>\n\n');

fprintf(h,'<center><h2>%s tag browser</h2></center>\n',experiment);
fprintf(h,'    <center>');
[last_ph,next_ph,last_tag,next_tag]=select_last_next(pt,ii,false);
if last_ph>0
  fprintf(h,'[<a href="../../%s/browser/%s.html">prev_phase</a>]',pt.tag{last_ph}(1:6),pt.tag{last_ph});
else
  fprintf(h,'[<b>prev_phase</b>]');
end

if last_tag>0
  fprintf(h,'[<a href="../../%s/browser/%s.html">prev_tag</a>]',pt.tag{last_tag}(1:6),pt.tag{last_tag});
else
  fprintf(h,'<b>[prev_tag]</b>');
end
fprintf(h,' [<a id="static_link" target="_parent" href="');
fprintf(h,'http://bicep0.caltech.edu/~bicep2/pipeline/reducplots/browser/?%s',pt.tag{ii});
fprintf(h,'">%s</a>] ',pt.tag{ii});
fprintf(h,'<SCRIPT LANGUAGE="JavaScript">\n');
fprintf(h,'<!--\n');
% fprintf(h,'document["static_link"].href="http://"+location.host+"/"+location.pathname+"?%s"\n',pt.tag{ii});
fprintf(h,'static_link_base="http://"+parent.location.host+parent.location.pathname+"?%s";\n',pt.tag{ii});
fprintf(h,'document.getElementById(''static_link'').href=static_link_base+"&"+parent.reduc_plot_type+"&"+parent.reduc_plot_rx+"&"+parent.reduc_plot_cut;\n');
fprintf(h,'//-->\n</SCRIPT>\n\n');
if next_tag<=length(pt.tag)
  fprintf(h,'[<a href="../../%s/browser/%s.html">next_tag</a>]',pt.tag{next_tag}(1:6),pt.tag{next_tag});
else
  fprintf(h,'<b>[next_tag]</b>');
end
if next_ph<=length(pt.tag)
  fprintf(h,'[<a href="../../%s/browser/%s.html">next_phase</a>]',pt.tag{next_ph}(1:6),pt.tag{next_ph});
else
  fprintf(h,'[<b>next_phase</b>]');
end
fprintf(h,'</center>\n');

fprintf(h,'<center>\n');
fprintf(h,'<table border="1" cellspacing="0" cellpadding="5">\n');
fprintf(h,'<tr><td align="center">\n');
fprintf(h,'  <b><a href="javascript:toggle_update_taginfo();">Tag info</a></b><br>\n');
fprintf(h,'<table border="0" cellspacing="0" cellpadding="5" id="taginfo" style="display:none;">\n');
fprintf(h,'<tr><td align="right"><b>Tag:</b></td><td align="left">%s</td>\n',pt.tag{ii});
fprintf(h,'    <td align="right"><b>Phase:</b></td><td align="left">%s</td>\n',pt.ph{ii});
fprintf(h,'    <td align="right"><b>Start:</b></td><td align="left">%s</td>\n',pt.tstart{ii});
fprintf(h,'<tr><td align="right"><b>Schedule:</b></td><td align="left">%s</td>\n',pt.sch{ii});
fprintf(h,'    <td align="right"><b>Scanset:</b></td><td align="left">%s</td>\n',pt.set{ii});
fprintf(h,'    <td align="right"><b>End:</b></td><td align="left">%s</td>\n',pt.tend{ii});
fprintf(h,'<tr><td align="right"><b>Type:</b></td><td align="left">%s</td>\n',pt.type{ii});
fprintf(h,'    <td align="right"><b>Elofs:</b></td><td align="left">%s</td>\n',pt.elofs{ii});
fprintf(h,'    <td align="right"><b>Good:</b></td><td align="left">%s</td>\n',pt.good{ii});
fprintf(h,'<tr><td align="right"><b>Comments:</b></td>\n');
fprintf(h,'    <td align="left" colspan="5">%s</td></tr>\n',pt.comment{ii});
fprintf(h,'</table></td></tr>\n');
fprintf(h,'</table>\n');
fprintf(h,'</center>\n');
fprintf(h,'<SCRIPT LANGUAGE="JavaScript">\n');
fprintf(h,'<!--\n');
fprintf(h,'document.getElementById(''taginfo'').style.display=parent.taginfo_disp;\n');
fprintf(h,'function toggle_update_taginfo(){\n');
fprintf(h,'  parent.toggle_taginfo();\n');
fprintf(h,'  document.getElementById(''taginfo'').style.display=parent.taginfo_disp;\n');
fprintf(h,'}\n');
fprintf(h,'//-->\n</SCRIPT>\n\n');

fprintf(h,'<center>\n');
fprintf(h,'<table border="1" cellspacing="0" cellpadding="5">\n');
fprintf(h,'<tr><td align="center">\n');
fprintf(h,'  <b><a href="javascript:toggle_update_tagstats();">Scanset stats</a></b><br>\n');
fprintf(h,'<table border="0" cellspacing="0" cellpadding="5" id="tagstats" style="display:none;">\n');
maxcol=6;
numfields=length(ks.fields)-1;
if numfields<=maxcol
  ncols=numfields;
  nrows=1;
else
  ncols=maxcol;
  nrows=ceil(numfields/maxcol);
end
for irow=1:nrows
  fprintf(h,'<tr>\n');
  for icol=1:ncols
    iField=1+(irow-1)+(icol-1)*nrows;
    if iField>numfields
      fprintf(h,'    <td>&nbsp;</td>\n');
      fprintf(h,'    <td>&nbsp;</td>\n');
      continue
    end
    fprintf(h,'    <td><b>%s:</b></td>\n',ks.fields{iField+1});
    if isempty(jj) || isnan(jj)
      tmp='';
    elseif iscell(ps.(ks.fields{iField+1}))
      tmp=ps.(ks.fields{iField+1}){jj};
    else
      tmp=ps.(ks.fields{iField+1})(jj);
    end
    if ischar(tmp)
      fprintf(h,'    <td>%s</td>\n',tmp);
    else
      fprintf(h,'    <td>%f</td>\n',tmp);
    end
  end
  fprintf(h,'</tr>\n');
end
fprintf(h,'</table></td></tr>\n');
fprintf(h,'</table>\n\n');
fprintf(h,'<SCRIPT LANGUAGE="JavaScript">\n');
fprintf(h,'<!--\n');
fprintf(h,'document.getElementById(''tagstats'').style.display=parent.tagstats_disp;\n');
fprintf(h,'function toggle_update_tagstats(){\n');
fprintf(h,'  parent.toggle_tagstats();\n');
fprintf(h,'  document.getElementById(''tagstats'').style.display=parent.tagstats_disp;\n');
fprintf(h,'}\n');
fprintf(h,'function set_reduc_plot_type(tt){\n');
fprintf(h,'  parent.set_reduc_plot_type(tt);\n');
fprintf(h,'  document.getElementById(''static_link'').href=static_link_base+"&"+parent.reduc_plot_type+"&"+parent.reduc_plot_rx+"&"+parent.reduc_plot_cut;\n');
fprintf(h,'}\n');
fprintf(h,'function set_reduc_plot_rx(rr){\n');
fprintf(h,'  parent.set_reduc_plot_rx(rr);\n');
fprintf(h,'  document.getElementById(''static_link'').href=static_link_base+"&"+parent.reduc_plot_type+"&"+parent.reduc_plot_rx+"&"+parent.reduc_plot_cut;\n');
fprintf(h,'}\n');
fprintf(h,'function set_reduc_plot_cut(cc){\n');
fprintf(h,'  parent.set_reduc_plot_cut(cc);\n');
fprintf(h,'  document.getElementById(''static_link'').href=static_link_base+"&"+parent.reduc_plot_type+"&"+parent.reduc_plot_rx+"&"+parent.reduc_plot_cut;\n');
fprintf(h,'}\n');

% Key bindings for switching plots
fprintf(h, '\n');
fprintf(h, '// Figure out if we are B2/Keck or later\n');
fprintf(h, 'var url = document.location.toString();\n');
fprintf(h, 'var isb2 = url.indexOf("bicep2") != -1;\n');
fprintf(h, 'var iskeck = url.indexOf("keck") != -1;\n');
fprintf(h, 'var rxmce = (isb2 || iskeck) ? "_rx" : "_mce";\n');
fprintf(h, '\n');
fprintf(h, '// Handle key presses\n');
fprintf(h, 'function plupdate_by_key(e) {\n');
fprintf(h, '  var key = undefined;\n');
fprintf(h, '  var handled = false;\n');
fprintf(h, '\n');
fprintf(h, '  if (e.key !== undefined) {\n');
fprintf(h, '    key = e.key;\n');
fprintf(h, '  } else if (e.keyIdentifier !== undefined) {\n');
fprintf(h, '    key = e.keyIdentifier; // Untested\n');
fprintf(h, '  } else if (e.keyCode !== undefined) {\n');
fprintf(h, '    key = e.keyCode; // Untested\n');
fprintf(h, '  }\n');
fprintf(h, '\n');
fprintf(h, '  if ("0" <= key && key <= "9") {\n');
fprintf(h, '    var rx = rxmce + (parseInt(key)-1).toString();\n');
fprintf(h, '    //console.log("Setting plot rx to " + rx);\n');
fprintf(h, '    set_reduc_plot_rx(rx);\n');
fprintf(h, '    handled = true;\n');
fprintf(h, '  }\n');
fprintf(h, '  else if (key == " ") {\n');
fprintf(h, '    //console.log("Setting plot rx to all");\n');
fprintf(h, '    set_reduc_plot_rx("");\n');
fprintf(h, '    handled = true;\n');
fprintf(h, '  }\n');
fprintf(h, '  else if (key == "!") {\n');
fprintf(h, '    //console.log("Setting plot cut to round1");\n');
fprintf(h, '    set_reduc_plot_cut("_round1");\n');
fprintf(h, '    handled = true;\n');
fprintf(h, '  }\n');
fprintf(h, '  else if (key == "@") {\n');
fprintf(h, '    //console.log("Setting plot cut to round2");\n');
fprintf(h, '    set_reduc_plot_cut("_round2");\n');
fprintf(h, '    handled = true;\n');
fprintf(h, '  }\n');
fprintf(h, '  else if (key == "#") {\n');
fprintf(h, '    //console.log("Setting plot cut to noise");\n');
fprintf(h, '    set_reduc_plot_cut("");\n');
fprintf(h, '    handled = true;\n');
fprintf(h, '  }\n');
fprintf(h, '  else if (key == "$") {\n');
fprintf(h, '    //console.log("Setting plot cut to none");\n');
fprintf(h, '    set_reduc_plot_cut("_none");\n');
fprintf(h, '    handled = true;\n');
fprintf(h, '  }\n');
fprintf(h, '\n');
fprintf(h, '  if (handled) {\n');
fprintf(h, '    plupdate();\n');
fprintf(h, '    e.preventDefault();\n');
fprintf(h, '  }\n');
fprintf(h, '}\n');
fprintf(h, '\n');
fprintf(h, '// Register the event handler\n');
fprintf(h, 'document.addEventListener("keypress", plupdate_by_key, true);\n');
fprintf(h,'//-->\n</SCRIPT>\n\n');

% Set up parts of the filename given as inputs to set_reduc_plot_url
tagdate=pt.tag{ii}(1:8); %date
ph=pt.ph{ii}; %phase
ss=sprintf('%02d',str2num(pt.set{ii})); %scanset
dk=sprintf('dk%03d',str2num(pt.dk{ii})); %deck angle

%All receivers on same plot / coadded
fprintf(h,'<center>\n');

fprintf(h,'<table border="0" cellspacing="1" cellpadding="1">\n');

fprintf(h,'<tr><th>Reduc plot type: </th>');
fprintf(h,'<td><a href="javascript:set_reduc_plot_type(''scans'');">scans</a></td>\n');
fprintf(h,'<td><a href="javascript:set_reduc_plot_type(''stddev'');">std</a></td>\n');
fprintf(h,'<td><a href="javascript:set_reduc_plot_type(''elnod'');">elnod</a></td>\n');
fprintf(h,'<td><a href="javascript:set_reduc_plot_type(''PSD'');">psd</a></td>\n');
fprintf(h,'<td><a href="javascript:set_reduc_plot_type(''NET'');">net</a></td>\n');
fprintf(h,'<td><a href="javascript:set_reduc_plot_type(''cov'');">cov</a></td>\n');
fprintf(h,'<td><a href="javascript:set_reduc_plot_type(''cal'');">cals</a></td>\n');
fprintf(h,'<td><a href="javascript:set_reduc_plot_type(''snow'');">snow</a></td>\n');
fprintf(h,'<td><a href="javascript:set_reduc_plot_type(''light'');">lights</a></td>\n');
fprintf(h,'<td><a href="javascript:set_reduc_plot_type(''dark_sd'');">darks</a></td>\n');
fprintf(h,'</tr>\n');
fprintf(h,'</table>\n');

if(strmatch(experiment,'keck') | strmatch(experiment,'Keck'))
  [p,ind] = get_array_info(pt.tag{ii});
  freqs = unique(p.band(ind.la));

  fprintf(h,'<table border="0" cellspacing="1" cellpadding="1">\n');
  fprintf(h,'<tr><th>rx: </th>\n');
  fprintf(h,'<td><a href="javascript:set_reduc_plot_rx('''');">all</a></td>\n');
  for ff=1:length(freqs)
    fprintf(h,'<td><a href="javascript:set_reduc_plot_rx(''_all%03d'');">all%03d</a></td>\n',freqs(ff),freqs(ff));
  end
  fprintf(h,'<td><a href="javascript:set_reduc_plot_rx(''_rx0'');">rx0</a></td>\n');
  fprintf(h,'<td><a href="javascript:set_reduc_plot_rx(''_rx1'');">rx1</a></td>\n');
  fprintf(h,'<td><a href="javascript:set_reduc_plot_rx(''_rx2'');">rx2</a></td>\n');
  if datenum(tagdate,'yyyymmdd')>=2012
    fprintf(h,'<td><a href="javascript:set_reduc_plot_rx(''_rx3'');">rx3</a></td>\n');
    fprintf(h,'<td><a href="javascript:set_reduc_plot_rx(''_rx4'');">rx4</a></td>\n');
  end
  fprintf(h,'</tr>\n');
  fprintf(h,'</table>\n');
end

if(strmatch(experiment,'bicep3'))
  fprintf(h,'<table border="0" cellspacing="1" cellpadding="1">\n');
  fprintf(h,'<tr><th>mce: </th>\n');
  fprintf(h,'<td><a href="javascript:set_reduc_plot_rx('''');">all</a></td>\n');
  fprintf(h,'<td><a href="javascript:set_reduc_plot_rx(''_mce0'');">mce0</a></td>\n');
  fprintf(h,'<td><a href="javascript:set_reduc_plot_rx(''_mce1'');">mce1</a></td>\n');
  fprintf(h,'<td><a href="javascript:set_reduc_plot_rx(''_mce2'');">mce2</a></td>\n');
  if datenum(tagdate,'yyyymmdd')>=datenum('20160101','yyyymmdd')
    fprintf(h,'<td><a href="javascript:set_reduc_plot_rx(''_mce3'');">mce3</a></td>\n');
  end
  fprintf(h,'</tr>\n');
  fprintf(h,'</table>\n');
end

fprintf(h,'<table border="0" cellspacing="1" cellpadding="1">\n');
fprintf(h,'<tr><th>cuts: </th>\n');
fprintf(h,'<td><a href="javascript:set_reduc_plot_cut('''');">noise</a></td>\n');
fprintf(h,'<td><a href="javascript:set_reduc_plot_cut(''_none'');">none</a></td>\n');
fprintf(h,'<td><a href="javascript:set_reduc_plot_cut(''_round1'');">round1</a></td>\n');
fprintf(h,'<td><a href="javascript:set_reduc_plot_cut(''_round2'');">round2</a></td>\n');
fprintf(h,'</tr>\n');
fprintf(h,'</table>\n');

fprintf(h,'<img src="%s" name="reduc">\n', ...
  fullfile('..',[pt.tag{ii} '_001.png']));
fprintf(h,'<SCRIPT LANGUAGE="JavaScript">\n');
fprintf(h,'<!--\n');
fprintf(h,'parent.set_reduc_plot_tag(''%s'',''%s'',''%s'',''%s'');\n', ...
tagdate,ph,ss,dk);
fprintf(h,'//-->\n</SCRIPT>\n');
fprintf(h,'</center>\n\n');

% Cuts
fprintf(h,'<center>Cuts: \n');
fprintf(h,'<a href="../cutmessages/%s_round1.txt">round1</a> \n',pt.tag{ii});
fprintf(h,'<a href="../cutmessages/%s_round2.txt">round2</a>\n',pt.tag{ii});
fprintf(h,'</center><p>\n\n');

% Box for plot help
fprintf(h,'<center>\n');
fprintf(h,'<table border="1" cellpadding="5" cellspacing="0" style="width:0; height:0;" id="plothelp_outer">\n');
fprintf(h,'<tr><td align="center">\n');
fprintf(h,'  <b><a href="javascript:toggle_update_plothelp();">Plot&nbsp;help</a></b><br>\n');
fprintf(h,'<table border="0" cellspacing="0" cellpadding="5" id="plothelp_inner" style="width:98%%; height:98%%; display:none;">\n');
fprintf(h,'<tr><td>\n');
fprintf(h,'<iframe width="98%%" height="95%%" src="" frameborder="0" id="plothelp_obj">\n');
fprintf(h,'</iframe></td></tr>\n');
fprintf(h,'</table>\n');
fprintf(h,'</td></tr>');
fprintf(h,'</table>\n');
fprintf(h,'</center>\n');

fprintf(h,'<SCRIPT LANGUAGE="JavaScript">\n');
fprintf(h,'<!--\n');
fprintf(h,'document.getElementById(''plothelp_inner'').style.display=parent.plothelp_disp;\n');
fprintf(h,'if (parent.plothelp_disp!="none"){\n');
fprintf(h,'  document.getElementById(''plothelp_outer'').style.width="80%%";\n');
fprintf(h,'  document.getElementById(''plothelp_outer'').style.height="75%%";\n');
fprintf(h,'  document.getElementById(''plothelp_obj'').src="../../browser/help_"+parent.reduc_plot_type+".html";\n');
fprintf(h,'};\n');
fprintf(h,'function toggle_update_plothelp(){\n');
fprintf(h,'  parent.toggle_plothelp();\n');
fprintf(h,'  document.getElementById(''plothelp_inner'').style.display=parent.plothelp_disp;\n');
fprintf(h,'if (parent.plothelp_disp=="none"){\n');
fprintf(h,'    document.getElementById(''plothelp_outer'').style.width="0";\n');
fprintf(h,'    document.getElementById(''plothelp_outer'').style.height="0";\n');
fprintf(h,'    document.getElementById(''plothelp_obj'').src="";\n');
fprintf(h,'  } else {\n');
fprintf(h,'    document.getElementById(''plothelp_outer'').style.width="80%%";\n');
fprintf(h,'    document.getElementById(''plothelp_outer'').style.height="75%%";\n');
fprintf(h,'    document.getElementById(''plothelp_obj'').src="../../browser/help_"+parent.reduc_plot_type+".html";\n');
fprintf(h,'  }\n');
fprintf(h,'}\n');
fprintf(h,'//-->\n</SCRIPT>\n\n');

fprintf(h,'</body>\n\n');
fprintf(h,'</html>');
fclose(h);

setpermissions(htmlfname);

return


function [last_ph next_ph last_tag next_tag]=select_last_next(pt,on_idx,do_show_all)
last_tag=0;
next_tag=length(pt.tag)+1;
last_ph=0;
next_ph=length(pt.tag)+1;

idx=on_idx;
while idx>0
  idx=idx-1;
  if idx<=0
    break
  end
  if ~do_show_all && strcmp(pt.ph{idx},'A')
    continue
  elseif ~do_show_all && strcmp(pt.set{idx},'0')
    continue
  end
  if last_tag==0
    last_tag=idx;
  end
  if strcmp(pt.ph{idx},pt.ph{on_idx}) && strncmp(pt.tag{idx},pt.tag{on_idx},8)
    continue
  end
  last_ph=idx;
  break
end

idx=on_idx;
while idx<=length(pt.tag)
  idx=idx+1;
  if idx>length(pt.tag)
    break
  end
  if ~do_show_all && strcmp(pt.ph{idx},'A')
    continue
  elseif ~do_show_all && strcmp(pt.set{idx},'0')
    continue
  end
  if next_tag>length(pt.tag)
    next_tag=idx;
  end
  if strcmp(pt.ph{idx},pt.ph{on_idx}) && strncmp(pt.tag{idx},pt.tag{on_idx},8)
    continue
  end
  next_ph=idx;
  break
end

return
