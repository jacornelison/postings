function [ra,dec]=arcvar_to_radec(d,racen,deccen)
% [ra,dec]=arcvar_to_radec(d,racen,deccen)
%
% Using data recorded in archive calculate the J2000 ra,dec
% pointing direction of the telescope at 100Hz rate
%
% Annoyingly the J2000 ra,dec of the current source is not
% recorded in the data stream. (The precession calc is done
% in navigator.c before sending to the rtc.)

if(~exist('raccen','var'))
  [racen,deccen]=src_coords(d.tracker.source);
end

% horiz_mount does not include encoder zero points but lockin.xxPos
% does - add them
d.tracker.horiz_mount=d.tracker.horiz_mount+d.tracker.encoder_off;

% horiz_mount does include scan and user offsets - remove them
d.tracker.horiz_mount=d.tracker.horiz_mount-d.tracker.scan_off;
d.tracker.horiz_mount=d.tracker.horiz_mount-d.tracker.horiz_off;

% Find pointing model correction (including refraction etc)
pointmod_az=d.tracker.horiz_mount(:,1)-d.tracker.horiz_geoc(:,1);
pointmod_el=d.tracker.horiz_mount(:,2)-d.tracker.horiz_geoc(:,2);
samprate=length(d.lockin.azPos)/size(d.tracker.horiz_mount,1);
pointmod_az =repmat(pointmod_az',samprate,1);  pointmod_az =cvec(pointmod_az);
pointmod_el =repmat(pointmod_el',samprate,1);  pointmod_el =cvec(pointmod_el);

% Remove pointing correction
az_geoc=d.lockin.azPos-pointmod_az;
el_geoc=d.lockin.elPos-pointmod_el;

% Find equat to horiz correction
% (this is good approx for locations close to Pole)
horiz2equat_az=d.tracker.horiz_geoc(:,1)-d.tracker.equat_geoc(:,1)*15;
horiz2equat_el=-d.tracker.horiz_geoc(:,2)-d.tracker.equat_geoc(:,2);
horiz2equat_az=repmat(horiz2equat_az',samprate,1);  horiz2equat_az =cvec(horiz2equat_az);
horiz2equat_el=repmat(horiz2equat_el',samprate,1);  horiz2equat_el =cvec(horiz2equat_el);

% Apply
ra =az_geoc-horiz2equat_az;
dec=-el_geoc-horiz2equat_el;

% Find precession (ra/dec) correction
precemod_ra=d.tracker.equat_geoc(:,1)*15-racen;
precemod_dec=d.tracker.equat_geoc(:,2)-deccen;
precemod_ra =repmat(precemod_ra',samprate,1);  precemod_ra =cvec(precemod_ra);
precemod_dec=repmat(precemod_dec',samprate,1); precemod_dec=cvec(precemod_dec);

% Apply to get back to J2000 coords
ra=ra-precemod_ra;
dec=dec-precemod_dec;

% Put ra into range 0-360
ind=ra<0; ra(ind)=ra(ind)+360;
ind=ra>360; ra(ind)=ra(ind)-360;

return
