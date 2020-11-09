function mkpng(filename,crop)
% function mkpng(filename,crop)
%
% Make a png image of the current figure

% Note that some extra complexity is required to make this work right,
% as 'getframe' doesn't seem to work on Odyssey...

if(~exist('crop','var'))
  crop=false;
end

% Capture screen figure
% x=getframe(fignum);

% Improved ugly hack follows.  Based on lengthy
% tinkering and tweaked with information from
% http://www.mathworks.com/support/solutions/en/data/1-16WME/?solution=1-16WME.


fignum=gcf;

% Horror: When you resize one figure, Matlab sometimes
% resizes ALL THE OTHERS!  Prevent this!
other_fignum=get(0,'Children');
for i=1:length(other_fignum)
  other_figpos{i}=get(other_fignum(i),'Position');
end

% Manipulate figure properties to cajole Matlab into saving
% with the right pixel count.
pp0=get(fignum,'PaperPosition');
ps0=get(fignum,'PaperSize');
pu0=get(fignum,'PaperUnits');
fu0=get(fignum,'Units');
vv0=get(fignum,'Visible');
pi0=get(0,'ScreenPixelsPerInch');

set(fignum,'Units','pixels');
fp0=get(fignum,'Position');
% In R2009a (and/or previous), manual fixups were required to make headless
% and headed figures appear [nearly] the same. With R2015a, most of these
% inconsistencies seem to have been resolved (probably corresponding to the
% major graphics subsystem rewrite in R2014b). We can distinguish these two
% cases simply by noting whether the graphics handle is a number (R2009-era)
% or not (an object handle after R2014b).
is_headless=strcmpi(get(fignum,'XDisplay'),'nodisplay') & isnumeric(fignum);
if is_headless
  % In headless mode, the -r100 option to print is
  % ignored.  Use the current pix-per-inch.
  pp1=fp0/pi0;
  % Sometimes in headless mode, corner can have neg.
  % coordinates on screen.  This breaks for paper pos.
  pp1(1)=0.5;
  pp1(2)=0.5;
else
  pp1=fp0/100;
end
% Make sure the hypothetical paper is plenty big enough.
ps1=[2*pp1(1)+pp1(3)+1 2*pp1(2)+pp1(4)+1];
set(fignum,'Visible','off');
set(fignum,'PaperSize',ps1);
set(fignum,'PaperUnits','inches','PaperPosition',pp1);

% If in headless mode, rescale all fonts so they'll
% look the same as in desktop mode
if is_headless
  of0=rescale_fonts(fignum,100/pi0);
end

% Save, cropping if necessary.
if(~crop)
  if is_headless
    print(fignum,'-dpng',filename);
  else
    print(fignum,'-dpng',filename,'-r100');
  end
else
  [fp fn fe]=fileparts(filename);
  if isempty(fe)
    fe='.png';
    filename=[filename fe];
  end
  tmpfname=fullfile(fp,[fn '.tmp' fe]);
  if is_headless
    print(fignum,'-dpng',tmpfname);
  else
    print(fignum,'-dpng',tmpfname,'-r100');
  end
  unix(sprintf('convert -trim -bordercolor white -border 20x20 %s %s',tmpfname,filename));
  system(['rm -f ' tmpfname]);
end

% Restore original figure properties.
set(fignum,'Units',fu0);
set(fignum,'PaperUnits',pu0);
set(fignum,'PaperPosition',pp0);
set(fignum,'PaperSize',ps0);
set(fignum,'Visible',vv0);
if is_headless
  unrescale_fonts(of0);
end

% Restore other figures' sizes.
for i=1:length(other_fignum)
  set(other_fignum(i),'Position',other_figpos{i});
end

return


% These functions are perhaps the ugliest of all ugly features in mkpng.
% When printing figures, Matlab scales points to pixels according to the
% current pixels-per-inch.  This is unavoidably different in desktop mode
% and headless mode; in fact, we can't effectively control it.  Thus, the
% only way to get the fonts to look the same is to adjust them by hand
% using the ratio of 100 dpi to the current pixels-per-inch number.  Must
% also remember to set them back at the end to the way they were at first.
function of=rescale_fonts(fignum,ratio)

of={};
hprop=get(fignum);
if isfield(hprop,'FontSize') && isfield(hprop,'FontUnits') && strcmp(hprop.FontUnits,'points')
  of={fignum, hprop.FontSize};
  if 1
    set(fignum,'FontSize',hprop.FontSize*ratio);
  else
    set(fignum,'FontUnits','pixels');
    set(fignum,'FontSize',hprop.FontSize*100/72);
  end
  % disp(['Rescaling fonts in ' hprop.Type ' ' num2str(fignum) ' from ' num2str(hprop.FontSize) ' to ' num2str(hprop.FontSize*ratio)]);
end
% Don't mess with children of a legend -- they're magically linked
% to the parent legend object.
if isfield(hprop,'Tag') && strcmp(hprop.Tag,'legend')
  return
end
child_fields={'Children','Title','XLabel','YLabel'};
for k=1:length(child_fields)
  if isfield(hprop,child_fields{k})
    for j=1:length(hprop.(child_fields{k}))
      of=[of; rescale_fonts(hprop.(child_fields{k})(j),ratio)];
    end
  end
end

return

% Put it back the way it was!
function unrescale_fonts(of)

for i=1:size(of,1)
  set(of{i,1},'FontUnits','points');
  set(of{i,1},'FontSize',of{i,2});
end

return

