function [x_B, y_B] = rtheta_to_xbyb(r, theta)
%
%Converts between (r, theta) and (x_B, y_B) beam map coordinates
%See http://bmode.caltech.edu/~spuder/analysis_logbook/analysis/20131005_sidelobe_coordinates/sidelobe_coordinates.pdf for coordinate definitions
%
%r, theta - in degrees
%
%rtod = 180 / pi;
%dtor = pi / 180;
d_B = 2 .* sind(r  ./ 2);
x_B = d_B .* cosd(theta );
y_B = d_B .* sind(theta );

end
