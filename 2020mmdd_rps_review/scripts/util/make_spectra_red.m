function make_spectra_red()
%makes a red spectrum for eb seperation using covariance matrix.

%get a template spectrum
spec_in='input_maps/official_cl/camb_planck2013_r0.fits';
spec=fitsread(spec_in, 'table');
info=fitsinfo(spec_in);
mod=load_cmbfast(spec_in);
hdr=info.PrimaryData.Keywords; xhdr=info.AsciiTable.Keywords; xhdr = check_units(xhdr);

%make an inverse l^2 spectrum 
inv_lsq=1./mod.l.^2;

inv_lsq(1)=0; %remove the inf for monopole

%now fill back in and save for BB=0 case
spec_out=spec;
spec_out{1}=inv_lsq; % TT
spec_out{2}=inv_lsq; % EE
spec_out{3}=zeros(size(spec{3})); %leave BB as zero
spec_out{4}=inv_lsq; % TE. is this ok for TE?


fitswrite_table(spec_out,hdr, xhdr, 'red_spectrum_EnoB.fits')

%now fill back in and save for EE=0 case
spec_out=spec;
spec_out{1}=inv_lsq; % TT
spec_out{2}=zeros(size(spec{2})); %leave EE as zero
spec_out{3}=inv_lsq; % BB
spec_out{4}=zeros(size(spec{4})); % TE. Since there is no E, TE should be zero as well.

% For BK analysis, we have applied the purification matrix
%  (composed of eigen mode solution separating pure E and pure B modes)
%  by using only the polarization-polarization components of pixel-pixel covariance matrices.
%  (see matrix/reduc_covmattheory.m for the calculation.)
% The TE power spectrum is used for TQ, TU, QT, UT blocks of the covariance matrices.
% For the analysis with constrained T, we consider only QQ, QU, UQ, UU blocks of
%  the covariance matrices and the reobservation matrices. Then the red TE spectrum doesn't do anything. 

fitswrite_table(spec_out,hdr, xhdr, 'red_spectrum_BnoE.fits')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = check_units(hdr)
for ii =1:size(hdr,1)
  if any(strfind(hdr{ii,1},'TUNIT'))
    if any(strfind(hdr{ii,2},'unknown')) hdr{ii,2} = 'K^2'; end
  end
end
return
