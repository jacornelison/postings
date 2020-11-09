function make_bmhtml(bmopt)
% make_bmhtml(bmopt)
% 
% Auto-generate an html file for beam maps
%
% INPUTS (should all be passed in with bmopt)
%
% expt:      'b2','keck','bicep3'
% rxNum:     For Keck: rx number, array of rx, or 'all' (default)
% filename:  string with timestamp of map file (i.e. '20150101_000000' or 
%            p = get_bm_info(year); filename = p.filename{bmnum};)
% plotdir:   directory in which images are saved, and in which index.html is
%            saved, default plots/filename
% pager:     0 (default) for traditional
%            1 for nifty pager stuff
% 
% STUFF FOR TOP OF HTML (nice but not necessary)
% author:    String to place on HTML as author (default ''), optional
% bmnum:     Number label to make bookkeeping easier (look at .csv)
% sched:     Name of schedule used to take map

% Parse bmopt and set defaults
expt = bmopt.expt;
if ~isfield(bmopt,'rxNum')
  switch expt
    case 'keck'
      rxNum = 0:4;
  end
else
  if strcmp(bmopt.rxNum,'all')
    rxNum = 0:4;
  else
    rxNum = bmopt.rxNum;
  end
end
filename = bmopt.filename;
if ~isfield(bmopt,'plotdir')
  plotdir = ['plots/' filename];
else
  plotdir = bmopt.plotdir;
end
if ~isfield(bmopt,'author')
  bmopt.author = '';
end
if ~isfield(bmopt,'bmnum')
  bmopt.bmnum = '';
end
if ~isfield(bmopt,'sched')
  bmopt.sched = '';
end
if ~isfield(bmopt,'pager')
  pager = 0;
else
  pager = bmopt.pager;
end

p = get_bm_info(str2num(filename(1:4)));

switch expt
  case {'b2','keck'}
    tilenum = [1 2 ; ...
               4 3];
    tilerot = [0 0 ; ...
               180 180];  
  case 'bicep3'
    switch bmopt.run
      case 'b3r6'
	tilenum = [ 5  6  7; ...
	            10 11 12; ...
		    15 16 17];
	tilerot = [ 0   180 180; ...
	            180 0   0; ...
		    180 180 0];
      case {'b3r8','b3r9'}
	tilenum = [ 0  14 9  4 -1; ...
	            19 15 10 5 1; ...
		    20 16 11 6 2; ...
		    0  17 12 7 3; ...
		    0  18 13 8 0];
	tilerot = [ 0   90  90  90  0; ...
	            90  90  270 270 90; ...
		    270 270 270 90  90; ...
		    0   0   180 270 270; ...
		    0   0   180 180 0];
    end
end

% Run the main function, once for each folder
switch expt
  case {'b2','bicep3'}
    gen_html(bmopt,plotdir,tilenum,tilerot,0,pager)
  case 'keck'
    for ii = 1:length(rxNum)
      plotdirrx = [plotdir '_rx' num2str(rxNum(ii))];
      gen_html(bmopt,plotdirrx,tilenum,tilerot,rxNum(ii),pager)
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gen_html(bmopt,plotdir,tilenum,tilerot,rxNum,pager)

mkdir(plotdir);
fid = fopen([plotdir '/index.html'],'wt');
fprintf(fid,'<!DOCTYPE html>\n');
fprintf(fid,'\n');
fprintf(fid,'<html lang="en">\n');
fprintf(fid,'  <head>\n');
shortdate = [bmopt.filename(1:4) '-' bmopt.filename(5:6)];
switch bmopt.expt
  case 'b2'
    fprintf(fid,['    <title>BICEP2 Beam Maps: ' shortdate '</title>\n']);
  case 'keck'
    fprintf(fid,['    <title>Keck Beam Maps: ' shortdate '</title>\n']);
  case 'bicep3'
    fprintf(fid,['    <title>BICEP3 Beam Maps: ' shortdate '</title>\n']);
end
fprintf(fid,'    <meta charset="utf-8" />\n');
% Add javascript pager stuff to top if required
make_pager_header(fid,pager)
fprintf(fid,'  </head>');
fprintf(fid,'\n');
fprintf(fid,'  <body style="color: #000000; background: #c8cfef;">\n');
fprintf(fid,'\n');
fprintf(fid,'    <center>\n');
switch bmopt.expt
  case 'b2'
    fprintf(fid,['    <h1>BICEP2 Beam Maps: ' shortdate '</h1>\n']);
  case 'keck'
    fprintf(fid,['    <h1>Keck Beam Maps: ' shortdate '</h1>\n']);
  case 'bicep3'
    fprintf(fid,['    <h1>BICEP3 Beam Maps: ' shortdate '</h1>\n']);
end
fprintf(fid, '    </center>\n');
fprintf(fid, '\n');
fprintf(fid,['    ' bmopt.author '<br>\n']);
fprintf(fid,['    ' datestr(now,'yyyy-mm-dd') '<br>\n']);
fprintf(fid, '    <hr>\n');
fprintf(fid, '\n');
switch bmopt.expt
  case 'keck'
    fprintf(fid,['    Rx' num2str(rxNum) '<br>\n']);
end
fprintf(fid,['    Beam Map Run ' num2str(bmopt.bmnum) ': ' bmopt.sched '<br>\n']);
fprintf(fid,['    t1 = ' bmopt.t1 '<br>\n']);
fprintf(fid,['    t2 = ' bmopt.t2 '<br>\n']);
fprintf(fid,'    <hr>\n');
fprintf(fid,'\n');

for polstate = {'A','B','dif'}
  
  % Header for Single polarization/difference
  fprintf(fid,'    <center>\n');
  fprintf(fid,'    <h2> %s Polarization </h2>\n', polstate{1});
  fprintf(fid,'    </center>\n');
  fprintf(fid,'\n');

  detrowlim = 1:8;
  detcollim = 1:8;
  % If Keck, choose between 100/{150/220} GHz tile sizes
  switch bmopt.expt
    case 'keck'
      switch bmopt.run
	case {'k2014','k2015'} 
	  switch rxNum
	    case {0,2}
	      detrowlim = 1:6;
	      detcollim = 1:6;
	  end
      end	  
  end
  
  if pager
    fprintf(fid,'  <table align="center">\n');
    fprintf(fid,'    <tr>\n');
    fprintf(fid,'      <td>\n');
  end
  % Start table for FP layout for each pol
  fprintf(fid,'    <table border="1" cellspacing="0" cellpadding="0" align="center">\n');
  for tilerow = 1:size(tilenum,1)
    fprintf(fid,'      <tr>\n');
    for tilecolumn = 1:size(tilenum,2)
      fprintf(fid,'        <td>\n');
      % Now for each tile
      if tilenum(tilerow,tilecolumn) > 0
	fprintf(fid,'          <table>\n');
	fprintf(fid,'            <tr>\n');
	fprintf(fid,'              <td colspan="6">Tile %u</td>\n',tilenum(tilerow,tilecolumn));
	if tilerot(tilerow,tilecolumn) == 0
	  for detrow = detrowlim
	    fprintf(fid,'            <tr>\n');
	    for detcolumn = detcollim
	      fprintf(fid,'              <td><a href="full/det_row_%u_col_%u_tile_%u_%s.png"><img src="small/det_row_%u_col_%u_tile_%u_%s_sm.png" width=35 class="gallery_%s"></a></td>\n',detrow,detcolumn,tilenum(tilerow,tilecolumn),polstate{1},detrow,detcolumn,tilenum(tilerow,tilecolumn),polstate{1},polstate{1});
	    end
	    fprintf(fid,'            </tr>\n');
	  end
	end
	if tilerot(tilerow,tilecolumn) == 90
	  for detcolumn = detrowlim
	    fprintf(fid,'            <tr>\n');
	    for detrow = detcollim
	      rownum = max(detrowlim) + 1 - detrow;
	      fprintf(fid,'              <td><a href="full/det_row_%u_col_%u_tile_%u_%s.png"><img src="small/det_row_%u_col_%u_tile_%u_%s_sm.png" width=35 class="gallery_%s"></a></td>\n',rownum,detcolumn,tilenum(tilerow,tilecolumn),polstate{1},rownum,detcolumn,tilenum(tilerow,tilecolumn),polstate{1},polstate{1});
	    end
	    fprintf(fid,'            </tr>\n');
	  end
	end
	if tilerot(tilerow,tilecolumn) == 180
	  for detrow = detrowlim
	    rownum = max(detrowlim) + 1 - detrow;
	    fprintf(fid,'            <tr>\n');
	    for detcolumn = detcollim
	      colnum = max(detcollim) + 1 - detcolumn;
	      fprintf(fid,'              <td><a href="full/det_row_%u_col_%u_tile_%u_%s.png"><img src="small/det_row_%u_col_%u_tile_%u_%s_sm.png" width=35 class="gallery_%s"></a></td>\n',rownum,colnum,tilenum(tilerow,tilecolumn),polstate{1},rownum,colnum,tilenum(tilerow,tilecolumn),polstate{1},polstate{1});
	    end
	    fprintf(fid,'            </tr>\n');
	  end
	end
	if tilerot(tilerow,tilecolumn) == 270
	  for detcolumn = detrowlim
	    colnum = max(detcollim) + 1 - detcolumn;
	    fprintf( fid, '            <tr>\n');
	    for detrow = detcollim
	      fprintf(fid,'              <td><a href="full/det_row_%u_col_%u_tile_%u_%s.png"><img src="small/det_row_%u_col_%u_tile_%u_%s_sm.png" width=35 class="gallery_%s"></a></td>\n',detrow,colnum,tilenum(tilerow,tilecolumn),polstate{1},detrow,colnum,tilenum(tilerow,tilecolumn),polstate{1},polstate{1});
	    end
	    fprintf(fid,'            </tr>\n');
	  end
	end
	fprintf( fid, '          </table>\n');
      elseif tilenum(tilerow, tilecolumn) == -1 % Pulse tube position
	% Draw a circle with PT in the middle
	fprintf(fid,'          <canvas id="circlecanvas" width="250" height="250"></canvas>\n');
	fprintf(fid,'          <script>\n');
	fprintf(fid,'            var canvas = document.getElementById("circlecanvas");\n');
	fprintf(fid,'            var context = canvas.getContext("2d");\n');
	fprintf(fid,'            context.arc(125,125,125,125,0,2*Math.PI, false);\n');
	fprintf(fid,'            context.stroke()\n');
	fprintf(fid,'          </script>\n');
      end
      fprintf(fid,'        </td>\n');
    end
    fprintf(fid,'      </tr>\n');
  end
  fprintf(fid,'    </table>\n');
  
  % And the actual paging element!
  if pager
    fprintf(fid,'      </td>\n');
    fprintf(fid,'      <td>\n');
    fprintf(fid,'        <div id="photoGallery">\n');
    fprintf(fid,'          <img id="photoLarge_%s" src="full/det_row_1_col_2_tile_1_%s.png" width="1000">\n',polstate{1},polstate{1});
    fprintf(fid,'        </div>\n');
    fprintf(fid,'      </td>\n');
    fprintf(fid,'    </table>\n');
  end
    
  
end

% Plot something about dark squids....

fprintf(fid,'</body>\n');
fprintf(fid,'\n');
fprintf(fid,'</html>');

fclose(fid);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_pager_header(fid,pager)

if pager
  % Maybe want a local jquery file?
  fprintf(fid,'    <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.12.0/jquery.min.js"></script>\n');
  fprintf(fid,'    <script>\n');
  fprintf(fid,'      $(function() {\n');
  fprintf(fid,'        $("a:has(img.gallery_A)").click(function() {\n');
  fprintf(fid,'          var largePath_A = $(this).attr("href");\n');
  fprintf(fid,'          var largePath_B = largePath_A.replace("_A","_B");\n');
  fprintf(fid,'          var largePath_dif = largePath_A.replace("_A","_dif");\n');
  fprintf(fid,'          $("#photoLarge_A").attr({src: largePath_A });\n');
  fprintf(fid,'          $("#photoLarge_B").attr({src: largePath_B });\n');
  fprintf(fid,'          $("#photoLarge_dif").attr({src: largePath_dif });\n');
  fprintf(fid,'          return false;\n');
  fprintf(fid,'        });\n');
  fprintf(fid,'      });\n');
  fprintf(fid,'\n');
  fprintf(fid,'      $(function() {\n');
  fprintf(fid,'        $("a:has(img.gallery_B)").click(function() {\n');
  fprintf(fid,'          var largePath_B = $(this).attr("href");\n');
  fprintf(fid,'          var largePath_A = largePath_B.replace("_B","_A");\n');
  fprintf(fid,'          var largePath_dif = largePath_B.replace("_B","_dif");\n');
  fprintf(fid,'          $("#photoLarge_A").attr({src: largePath_A });\n');
  fprintf(fid,'          $("#photoLarge_B").attr({src: largePath_B });\n');
  fprintf(fid,'          $("#photoLarge_dif").attr({src: largePath_dif });\n');
  fprintf(fid,'          return false;\n');
  fprintf(fid,'        });\n');
  fprintf(fid,'      });\n');
  fprintf(fid,'\n');
  fprintf(fid,'      $(function() {\n');
  fprintf(fid,'        $("a:has(img.gallery_dif)").click(function() {\n');
  fprintf(fid,'          var largePath_dif = $(this).attr("href");\n');
  fprintf(fid,'          var largePath_A = largePath_dif.replace("_dif","_A");\n');
  fprintf(fid,'          var largePath_B = largePath_dif.replace("_dif","_B");\n');
  fprintf(fid,'          $("#photoLarge_A").attr({src: largePath_A });\n');
  fprintf(fid,'          $("#photoLarge_B").attr({src: largePath_B });\n');
  fprintf(fid,'          $("#photoLarge_dif").attr({src: largePath_dif });\n');
  fprintf(fid,'          return false;\n');
  fprintf(fid,'         });\n');
  fprintf(fid,'      });\n');
  fprintf(fid,'    </script>\n');
end

return