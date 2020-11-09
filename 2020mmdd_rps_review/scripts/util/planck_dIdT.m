function dIdT=planck_dIdT(nu,T)
% dIdT=planck_dIdT(nu,T)
%
% Partial derivative of the Planck function with respect
% to temperature at given temperature and frequency.
%
% Input in Hz and K
% Output in W m^-2 Hz^-1 ster^-1 K^-1
%
% eg: To reproduce the plot Carlstrom uses everywhere:
%
% nu=0:1e9:500e9; Tcmb=2.726; y=1e-4;
% fx=sz_fx(nu); delTsz=Tcmb.*fx*y;
% dIdT=planck_dIdT(nu,Tcmb);
% delIsz=delTsz.*dIdT; plot(nu,delIsz.*1e7*1e-4);

  [h,k,c]= constants() ;
  x=(h*nu)./(k*T);

dIdT=((2*h^2*nu.^4)./(c^2*k*T.^2)).*(exp(x)./(exp(x)-1).^2);
