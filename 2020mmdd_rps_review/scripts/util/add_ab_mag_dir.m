function add_ab_mag_dir(fname)
% add_ab_mag_dir('csv_filename.csv')
%
% Add fields AB_mag and AB_dir to a beams CSV file. These are calculated from the Rs
% and THETAs in the CSV file itself. AB_mag is the great circle separation in
% degrees. AB_dir is the azimuth of the great circle connecting A and B at A (and
% pointing towards B.) If the fields exist, they are replaced. The respective fields
% are identical for the A and B member of each pair.

% get the ab offset parameters
[p,k]=ParameterRead(fname);

% get the array indices
[dum,ind]=get_array_info('20110505'); % random good day for Keck and BICEP2

% RA/DEC of individual detectors
[LAT,LON]=reckon(0,0,p.r,p.theta);

% Distance and azimuth of great circle connecting pairs
[D,AZ]=distance(LAT(ind.a),LON(ind.a),LAT(ind.b),LON(ind.b));

ab_mag=zeros(size(p.r));
ab_dir=zeros(size(p.r));

% Each member of the pair gets the same ab_mag and ab_dir
ab_mag(ind.a)=D;
ab_mag(ind.b)=D;
ab_dir(ind.a)=AZ;
ab_dir(ind.b)=AZ;

% Get fieldnames already present
fn=fieldnames(p);

if any(strcmp(fn,'ab_mag'))
  p.ab_mag=ab_mag;
else
  p.ab_mag=ab_mag;
  k.fields=[k.fields,'ab_mag'];
  k.units=[k.units,'deg'];
  k.formats=[k.formats,'float'];
end

if any(strcmp(fn,'ab_dir'))
  p.ab_dir=ab_dir;
else
  p.ab_dir=ab_dir;
  k.fields=[k.fields,'ab_dir'];
  k.units=[k.units,'deg'];
  k.formats=[k.formats,'float'];
end

ParameterWrite(fname,p,k);

return


