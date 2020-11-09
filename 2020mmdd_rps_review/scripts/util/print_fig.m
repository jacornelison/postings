%function print_fig(name, fig_handle, high_res, nocomp)
%
% Creates pdf in vector format and high resolution png figure 
%
%  Inputs: 
%    
%   name: filename of output figures will be name.pdf and name.png
%   fig_handle: figure handle 
%   high_res  : =1 (default) png will have crazy high resolution
%               =0 lower resolution
%               =2 lowest resolution
%   nocomp    : either 0 (default) or 1. 
%               set to 1 will demand the eps to pdf step
%               does not compress images. This is useful for maps.
%               Warning: file names with wildcards might not work with
%               this option
%%%%%%%%%%%%%%%%%%%%%%%%
function print_fig(name, fig_handle, high_res, nocomp)

if(~exist('fig_handle', 'var')|| isempty(fig_handle))
  fig_handle=gcf;
end

if(~exist('high_res', 'var')|| isempty(high_res))
  high_res=1;
end

if ~exist('nocomp', 'var') || isempty(nocomp)
  nocomp=0;
end


switch high_res
  case 0
    density='200'; %resolution of png in pixels per inch
  case 1
    density='500'; %resolution of png in pixels per inch
  case 2
    density='80'; %lowest res
end

set(fig_handle,'PaperPositionMode','auto');  %set printing to window size
set(fig_handle, 'Renderer', 'painters')
set(fig_handle, 'RendererMode', 'manual')
print(fig_handle, '-depsc2', '-Painters',name); %need painters for vector
fix_lines([name '.eps'])
if ~nocomp
  system_safe(['epstopdf ' name '.eps']);
else
  %for some reason, gs command doesn't work with ~
  if regexp(name, '~/')
    h=[getenv('HOME') '/'];
    name=strrep(name, '~', h);
  end
   system_safe(['gs -q -dBATCH -dNOPAUSE -dSAFER -dEPSCrop -sDEVICE=pdfwrite -sOutputFile=' name '.pdf -dAutoFilterMonoImages=false -dAutoFilterGrayImages=false -dAutoFilterColorImages=false -sColorImageFilter=FlateEncode -sGrayImageFilter=FlateEncode -sCompressPages=false -dPreserveHalftoneInfo=true ' name '.eps']);
end
%system_safe(['rm ' name '.eps']);  %maybe we want to keep eps for 
                                    %map images, which look better if ...
				    %you convert to pdf using an ...
				    %external program
				      

system_safe(['convert -density ' density  ' -units pixelsperinch ' name '.pdf ' name '.png']);

return
