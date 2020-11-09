function mkgif(filename,crop)
% function mkgif(filename,crop)
%
% Make a gif image of the current figure

if(~exist('crop','var'))
  crop=false;
end

% Capture screen figure
x=getframe(gcf);

% write ppm file
imwrite(x.cdata,'temp.ppm');

% Compose conversion command
if(~crop)
  cmd=sprintf('ppmquant 256 temp.ppm | ppmtogif > %s',filename);
else
  cmd=sprintf('pnmcrop -margin=20 temp.ppm | ppmquant 256 | ppmtogif > %s',filename);
end

% excute the command
unix(cmd);

% remove the temporary file
system('rm -f temp.ppm*');

return
