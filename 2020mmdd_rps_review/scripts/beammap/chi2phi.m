function [ phi ] = chi2phi(r,theta,chi)

% Loop over detectors
phi = zeros(size(r));
for i=1:length(phi)
    ri = r(i);
    th = theta(i);
    chii = chi(i);
    phi(i) = atan2(cosd(ri)*sind(th)+cosd(th)*tand(chii),cosd(ri)*cosd(th)-sind(th)*tand(chii))*180/pi;
end