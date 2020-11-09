function map=interpft2(map,s)
% map=interpft2(map,s)
%
% upsample array map to size s using Matlab provided interpft.m
% function
%
% see test code below
%

% cannot take ft when there are nans
map(isnan(map))=0;

map=interpft(map,s(1),1);
map=interpft(map,s(2),2);

return


% get a small test grid
m=get_map_defn('fttest_even');
% fill it with random numbers
map.Q=randn(m.ny,m.nx); map.U=randn(m.ny,m.nx);
ad=calc_ad2([m.xdos,m.ydos],[m.nx,m.ny]);
% calc the FTs as normal
ft=calc_map_fts(map,ad,1,[2],'normal',false);
% if we use interpft interatively...
mapift.Q=interpft2(map.Q,[m.ny,m.nx]*3); mapift.U=interpft2(map.U,[m.ny,m.nx]*3);
% we get a map which has the superset property...
map.Q
mapift.Q(1:3:end,1:3:end)
% and this can show us how to split the Nyquist elements in pad_ft...
ad2=calc_ad2([m.xdos,m.ydos],[m.nx,m.ny]*3);
ftift=calc_map_fts(mapift,ad2,1,[2],'normal',false);
ft.Q
z=ftift.Q; z(abs(z)<0.00001)=0
% implementing these emprically derived rules in pad_ft...
ftp=pad_ft(ft,[m.ny,m.nx]*3);
mapp.Q=f2i(ad,ftp.Q);
% we show that the superset property is reproduced
map.Q
mapp.Q(1:3:end,1:3:end)
