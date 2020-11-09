function mapo=maps_make_map(map,i,cflag)
% mapo=maps_make_map(map,i,cflag)
%
% Convert from tot signal and nhits to plotable map

map.map=map.tot(:,:,i)./map.nhit(:,:,i);
mapo=maps_conv_map(map,cflag);

return
