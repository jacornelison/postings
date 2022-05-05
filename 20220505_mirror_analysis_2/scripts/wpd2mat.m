function wpd2mat()

fname = 'z:/dev/moon_analysis/wpd_data_gary_et_al.csv';
f = fopen(fname);
Tdata = textscan(f,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',','HeaderLines',2);
%Tdata = csvread(fname,1);
T = [280:-10:200 180:-20:120 110 190:-20:130];

[phase, lat, BT] = deal([]);

for ind = 1:(length(Tdata)/2)
        phase = [phase; Tdata{2*(ind-1)+1}];
        lat = [lat; Tdata{2*ind}];
        BT = [BT; ones(size(Tdata{1}))*T(ind)];
end

phase = wrapTo180(phase+180);

ind = ~isnan(phase) & ~isnan(lat);
phase = phase(ind);
lat = lat(ind);
BT = BT(ind);

% Expand to account for wraps in phase
phase = [phase; phase+360; phase-360];
lat = [lat; lat; lat];
BT = [BT; BT; BT];

% and the mirroring of lat
phase = [phase; phase];
lat = [lat; -lat];
BT = [BT; BT];

fig = figure(3);
clf;
scatter(phase,lat,14,BT,'filled')
colorbar()
colormap default
grid on

save('z:/dev/moon_analysis/moon_temp_data.mat','BT','lat','phase')


