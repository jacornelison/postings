function write_ascii_map(ft, fname, nx, ny, type)
% write_ascii_map(map, fname, purified, nx, ny, type)
% 
% Save map to ASCII file
%
% ft       : FT data to be saved
% fname    : output file name
% nx,ny    : number of x, y grids
% type     : 't'        -> write FT of T
%          : 'eb'       -> write E/B modes
%
% examples:
% write_ascii_map(ft,'test.dat',236,100,'eb');
%

% ---- select output type ----
switch type
case 't' % output FT of temperature

  % 2d to 1d array in fortran order
  T = reshape(ft.T,nx*ny,1);

  % nan to zeros
  T(isnan(T)) = 0;

  % save into file
  fID = fopen(fname,'w');
  for i = 1:nx*ny
    fprintf(fID,'%12.5E %12.5E\n',real(T(i)),imag(T(i)));
  end
  fclose(fID);


% output E/B modes
case 'eb' 

  % 2d to 1d array in fortran order
  E = reshape(ft.E,nx*ny,1);
  B = reshape(ft.B,nx*ny,1);

  % nan to zeros
  E(isnan(E)) = 0;
  B(isnan(B)) = 0;

  % save into file
  fID = fopen(fname,'w');
  for i = 1:nx*ny
    fprintf(fID,'%12.5E %12.5E %12.5E %12.5E\n',real(E(i)),imag(E(i)),real(B(i)),imag(B(i)));
  end
  fclose(fID);


% output FT of T and E/B modes
case 'teb' 

  % 2d to 1d array in fortran order
  T = reshape(ft.T,nx*ny,1);
  E = reshape(ft.E,nx*ny,1);
  B = reshape(ft.B,nx*ny,1);

  % nan to zeros
  T(isnan(T)) = 0;
  E(isnan(E)) = 0;
  B(isnan(B)) = 0;

  % save into file
  fID = fopen(fname,'w');
  for i = 1:nx*ny
    fprintf(fID,'%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n',real(T(i)),imag(T(i)),real(E(i)),imag(E(i)),real(B(i)),imag(B(i)));
  end
  fclose(fID);


% no output
otherwise
  disp('no output')

end

