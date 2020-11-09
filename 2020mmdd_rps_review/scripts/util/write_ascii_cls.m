function write_ascii_cls(almfile,clfile)
% write_ascii_cls(almfile,clfile)
% 
% Save Cls to ASCII file
%
% almfile : fits filename of alm to be used for computing Cls
% clfile  : output ascii filename
%
% examples:
% write_ascii_cls('test.fits','test.txt');
%

disp('read alms')
almT=read_fits_alms(almfile,1);
almE=read_fits_alms(almfile,2);
almB=read_fits_alms(almfile,3);
disp('alm -> cl')
[l,TT]=alm2cl(almT);
[l,EE]=alm2cl(almE);
[l,BB]=alm2cl(almB);
fID = fopen(strcat(clfile),'w');
for i = 1:length(l)
  fprintf(fID,'%12.5E %12.5E %12.5E %12.5E\n', l(i), TT(i), EE(i), BB(i));
end
fclose(fID);

return

