function write_ascii_map(map, fname, nx, ny, type, amp)
% write_ascii_map(map, fname, purified, nx, ny, type)
% 
% Save map to ASCII file
%
% map      : map data to be saved
% fname    : output file name
% nx,ny    : number of x, y grids
% type     : 't'        -> write T map
%          : 'qu'       -> write normal Q/U maps
%          : 'qu_purif' -> write purified Q/U maps
% amp      : overall amplitude of map (default=1). Currently only work for temperature
%
% examples:
% write_ascii_map(map,'test.dat',236,100,'qu_purif');
%

% overall amplitude of map
if(~exist('amp','var'))   amp = 1.0; end;


% ---- select output type ----
switch type
case 't' % output temperature map

  % 2d to 1d array in fortran order
  T = reshape(map.T,nx*ny,1);

  % amplify
  if (amp~=1.0) T = T*amp; end;

  % nan to zeros
  T(isnan(T)) = 0;

  % save into file
  fID = fopen(fname,'w');
  for i = 1:nx*ny
    fprintf(fID,'%12.5E\n', T(i));
  end
  fclose(fID);


% output normal Q/U maps
case 'qu' 

  % 2d to 1d array in fortran order
  Q = reshape(map.Q,nx*ny,1);
  U = reshape(map.U,nx*ny,1);

  % nan to zeros
  Q(isnan(Q)) = 0;
  U(isnan(U)) = 0;

  % save into file
  fID = fopen(fname,'w');
  for i = 1:nx*ny
    fprintf(fID,'%12.5E %12.5E\n', Q(i), U(i));
  end
  fclose(fID);


% output purified Q/U maps
case 'qu_proj' 

  % 2d to 1d array in fortran order
  QB = reshape(full(map.QprojB),nx*ny,1);
  UB = reshape(full(map.UprojB),nx*ny,1);
  QE = reshape(full(map.QprojE),nx*ny,1);
  UE = reshape(full(map.UprojE),nx*ny,1);

  % nan to zeros
  QB(isnan(QB)) = 0;
  UB(isnan(UB)) = 0;
  QE(isnan(QE)) = 0;
  UE(isnan(UE)) = 0;

  % save into file
  fID = fopen(fname,'w');
  for i = 1:nx*ny
    fprintf(fID,'%12.5E %12.5E %12.5E %12.5E\n', QE(i), UE(i), QB(i), UB(i));
  end
  fclose(fID);

case 'tqu' % output temperature, Q and U map

  % 2d to 1d array in fortran order
  T = reshape(map.T,nx*ny,1);
  Q = reshape(map.Q,nx*ny,1);
  U = reshape(map.U,nx*ny,1);

  % amplify
  if (amp~=1.0) T = T*amp; end;

  % nan to zeros
  T(isnan(T)) = 0;
  Q(isnan(Q)) = 0;
  U(isnan(U)) = 0;

  % save into file
  fID = fopen(fname,'w');
  for i = 1:nx*ny
    fprintf(fID,'%12.5E %12.5E %12.5E\n', T(i), Q(i), U(i));
  end
  fclose(fID);


% no output
otherwise
  disp('no output')

end

