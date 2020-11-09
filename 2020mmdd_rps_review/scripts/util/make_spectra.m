function make_spectra(dir, r)
% function make_spectra(dir, r)
% makes spectra for planck LCMD
% vary tensor to scal ratio, r
% this needs to be run in the camb folder
%  
% dir   : output directory for files
% r    : tensor to scalar ratio
%  
% make_spectra('/n/home01/jetolan/biceppipe/input_spectra/',1);
% make_spectra('/n/home01/jetolan/biceppipe/input_spectra/',0.1);
% make_spectra('/n/home01/jetolan/biceppipe/input_spectra/',0);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %what input parameters are used in run_camb?
  params='planck2013';
 
  if ~exist('dir','var')
  dir='';
  end

  if ~exist('r','var')
  r=[];
  end

  if(isempty(dir))
  dir='';
  end
  
  if(isempty(r))
  r=0.1
  end

%make a string for this r  
  strr=num2str(r);
  pp=strfind(strr, '.');
  if ~isempty(pp)
  strrp=strrep(strr, '.', 'p');
else
  strrp=strr
  end
  
  %factor to multiply EE by for 'noE' case
  fac=1d-20 
  

  
%makes spectra using camb
%------------------------------------------------------
%with E, no lensing
%------------------------------------------------------
%%%%%%%%

  %run camb
run_camb([dir 'camb_' params '_r' strrp],'T','T', strr) 

  %read in spec
spec=fitsread([dir 'camb_' params '_r' strrp '_scalcls.fits'], 'table');

%get hdr info from scalar file
info=fitsinfo([dir 'camb_' params '_r' strrp '_scalcls.fits']);
hdr=info.PrimaryData.Keywords; xhdr=info.AsciiTable.Keywords; xhdr = check_units(xhdr);

%write our newly created spectra:
fitswrite_table(spec, hdr, xhdr, [dir 'camb_' params '_r' strrp '.fits'])

plot_spectra([dir 'camb_' params '_r' strrp '.fits'])


%%%%%%%%%
%------------------------------------------------------
%with E, with lensing
%------------------------------------------------------

%read camb file created above that has lensing:
spec=dlmread([dir 'camb_' params '_r' strrp '_lensedtotcls.dat']);

%get hdr info from scalar file
info=fitsinfo([dir 'camb_' params '_r' strrp  '_scalcls.fits']);
hdr=info.PrimaryData.Keywords; xhdr=info.AsciiTable.Keywords; xhdr = check_units(xhdr);

%massage spec so that it can be substittuted into .fit header
for i=1:4
  os=zeros(2001,1);
  ds=spec(:,i+1);                     %grab each spectra
  ds=ds.*(2*pi)*1d-12;                %convert units to K
  ds=ds./(spec(:,1).*(spec(:,1)+1));  %divide by l(l+1) as this is standard format
  os(spec(:,1)+1)=ds;                 %shift by one since first entry is l=0
  fin{i}=os;
end

%write our newly created spectra:
fitswrite_table(fin, hdr, xhdr, [dir 'camb_' params '_r' strrp '_lensing.fits'])

plot_spectra([dir 'camb_' params '_r' strrp '_lensing.fits'])

%------------------------------------------------------
%no E, no lensing
%------------------------------------------------------
%Remove EE power

%read in normal power
spec_in=fitsread([dir 'camb_' params '_r' strrp '.fits'], 'table');

%get hdr info from scalar file
info=fitsinfo([dir 'camb_' params '_r' strrp '.fits']);
hdr=info.PrimaryData.Keywords; xhdr=info.AsciiTable.Keywords; xhdr = check_units(xhdr);

%combine the two
for i=1:4
  ss=spec_in{i};
  
  if i==2  %EE is second entry
    fin{i}=ss.*fac;
  elseif i==4   %TE is 4th
    fin{i}=ss.*sqrt(fac);
    else
  fin{i}=ss;
  end
  
end

%write our newly combined spectra:
fitswrite_table(fin, hdr, xhdr, [dir 'camb_' params '_r' strrp  '_noE.fits'])

plot_spectra([dir 'camb_' params '_r' strrp '_noE.fits'])

%------------------------------------------------------
%no E, with lensing
%------------------------------------------------------
%Remove EE power

%read in normal power
spec_in=fitsread([dir 'camb_' params '_r' strrp '_lensing.fits'], 'table');

%get hdr info from scalar file
info=fitsinfo([dir 'camb_' params '_r' strrp '.fits']);
hdr=info.PrimaryData.Keywords; xhdr=info.AsciiTable.Keywords; xhdr = check_units(xhdr);


%combine the two
for i=1:4
  ss=spec_in{i};
  
  if i==2  %EE is second entry
    fin{i}=ss.*fac;
  elseif i==4   %TE is 4th
    fin{i}=ss.*sqrt(fac);
    else
  fin{i}=ss;
  end
  
end

%write our newly combined spectra:
fitswrite_table(fin, hdr, xhdr, [dir 'camb_' params '_r' strrp  '_lensing_noE.fits'])

plot_spectra([dir 'camb_' params '_r' strrp '_lensing_noE.fits'])

%------------------------------------------------------
%r=1, no scalars
%------------------------------------------------------
%%%%%%%%
% run camb
run_camb([dir 'camb_' params '_r' strrp '_noscal'],'F','T', strr) 
  
% read in spec
spec=fitsread([dir 'camb_' params '_r' strrp '_noscal_scalcls.fits'], 'table');

%get hdr info from scalar file
info=fitsinfo([dir 'camb_' params '_r' strrp '_noscal_scalcls.fits']);
hdr=info.PrimaryData.Keywords; xhdr=info.AsciiTable.Keywords; xhdr = check_units(xhdr);

%write our newly created spectra:
fitswrite_table(spec, hdr, xhdr, [dir 'camb_' params '_r' strrp '_noscal.fits'])

plot_spectra([dir 'camb_' params '_r' strrp '_noscal.fits'])


%------------------------------------------------------
% spectra for scan for n_t done in reduc_rlim
%------------------------------------------------------
% the tensor pivot scale here shifted away from 0.002/Mpc to
% break the degeneracy between scaling and slope:
p_t=0.009
for n_t=linspace(-2.5,2.5,51)
  filename=[dir 'camb_planck2013_r0p1_nt' num2str(n_t) '_pt' num2str(p_t)];
  run_camb(filename,'T','T','0.1', num2str(n_t), num2str(p_t))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot


function plot_spectra(file)

spec=load_cmbfast(file);


setwinsize(gcf, 1000,1000)
subplot(3,1,1)
plot(spec.Cs_l(1:500,1))

title('TT')
ylabel('l(l+1)C_l/2\pi (\muK)')
xlabel('l')

subplot(3,1,2)
plot(spec.Cs_l(1:500,3))

title('EE')
ylabel('l(l+1)C_l/2\pi (\muK)')
xlabel('l')

subplot(3,1,3)
plot(spec.Cs_l(1:500,4))

title('BB')
ylabel('l(l+1)C_l/2\pi (\muK)')
xlabel('l')

print('-dpng', strtok(file, '.'))

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = check_units(hdr)
for ii =1:size(hdr,1)
  if any(strfind(hdr{ii,1},'TUNIT'))
    if any(strfind(hdr{ii,2},'unknown')) hdr{ii,2} = 'K^2'; end
  end
end
return

