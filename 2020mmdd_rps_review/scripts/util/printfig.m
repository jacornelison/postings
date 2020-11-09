function printfig(f, filename, type, res, bclr)

% PRINTFIG prints image of figure to file
%
% printfig(f, filename, type, backgroundcolor) prints image of figure 'f' to file 
% 'filename' of optional type 'type'. Default is 'gif', other options
% are 'ps', 'eps' for postscript, 'pp' for a large gif generated from 
% postscript suitable for powerpoint. See function imwrite for other 
% type options. Background of figure is set to  white for printing, 
% then reset to original color. Optional 4th argument specifies background
% color for printing.

if nargin < 3
  type = 'gif';
end

if nargin  > 3
  xsize = res;
else
  xsize= 1000;
end

colr = get(f,'Color');
set(f,'Color',colr);

if strcmp(type,'ps') | strcmp(type,'eps') | strcmp(type,'pp')
  % Preserve aspect ratio of figure, and rescale to fill page. 
%  set(f,'PaperUnits','inches');
  set(f,'PaperPositionMode','auto');
%  pappos = get(f,'PaperPosition');
  % Now rescale to fit on page, with boundary 6w by 9h, papersize 8.5 x 11
%  if pappos(4)/pappos(3) > 9/6
%    % Scale image.
%    pappos(3:4) = pappos(3:4)/pappos(4)*9;
%  else
%    pappos(3:4) = pappos(3:4)/pappos(4)*6;
%  end
%  % Set lower left corner offset.
%  pappos(1:2) = ([8.5 11] - pappos(3:4))/2;
%  set(f,'PaperPosition',pappos);  
  if strcmp(type,'ps')
    eval(['print -f',num2str(f),' -dpsc2 ', filename],'error(lasterr)');
    return;
  elseif strcmp(type,'eps')
    eval(['print -f',num2str(f),' -depsc2 ', filename],'error(lasterr)');
    return;
  elseif strcmp(type,'pp')
    tmpfile = tempname;
    tmpfile = [tmpfile,'.eps'];
    print(['-f',num2str(f)],'-depsc2',tmpfile)
    [err str] = unix(['pstopnm -portrait -xborder 0.05 -yborder 0.05 -xsize ',num2str(xsize), ' ',tmpfile]);
    if err
      error(str);
    end
    % [err str] = unix(['ppmtogif ', tmpfile(6:end),'001.ppm > ',filename]);
    [err str] = unix(['ppmtogif ', tmpfile,'001.ppm > ',filename]);
    if err
      error(str);
    end
    % [err str] = unix(['\rm ', tmpfile,' ',tmpfile(6:end),'001.ppm']);
    [err str] = unix(['\rm ', tmpfile,' ',tmpfile,'001.ppm']);
    if err
      error(str);
    end
    
    set(f,'PaperPositionMode','manual')
    return;
  end
end


if strcmp(type,'png')
  % this re-renders the figure, but at least preserves the same size (in pixels) as
  % on the screen
  set(figure(f),'paperpositionmode','auto')
  eval(['print -f',num2str(f),' -dpng -r90 ', filename],'error(lasterr)');
  return;
end  
  
% bring figure to foreground
figure(f);
drawnow

% Get background color and set to white
colr = get(f,'Color');
if nargin == 5
  set(f,'Color',bclr);
else
  set(f,'Color','white');
end

img = getframe(f);

if strcmp(type,'gif')
  tmpfile = tempname;
  imwrite(img.cdata,tmpfile,'tif');
  [err str] = unix(['tifftopnm ', tmpfile, ...
	'| ppmquant 256 | ppmtogif >! ', filename]);
  if err
    error(str);
  end
  [err str] = unix(['\rm ', tmpfile]);
  if err
    error(str);
  end
else  
  imwrite(img.cdata,filename,type);
end

set(f,'Color',colr);
