function mapo=maps_conv_map(map,cflag)
% mapo=maps_conv_map(map,cflag)
%
% convolve map with gaussian for plotting

w=0.04; % sigma of gaussian convolution

if(cflag==1)
  delx=(map.x_tic(2)-map.x_tic(1))*cos(50*pi/180);
  dely=map.y_tic(2)-map.y_tic(1);
  x=-3*w:delx:3*w;
  y=-3*w:dely:3*w;
  [xg,yg]=meshgrid(x,y);
  beam=gauss2([1,0,0,w],xg,yg);
  beam=beam./sum(beam(:));
  mapo=conv2(map.map,beam,'same');
else
  mapo=map.map;
end

return
