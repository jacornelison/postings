function d=arcvar_to_azel(d)
% d=arcvar_to_azel(d)
%
% Using data recorded in archive calculate the az/el and offsets
% BICEP CMB fields are "CENTER" sources in source.cat which are
% specified as RA/Dec - a "track" command is used on these but
% horiz_topo does not appear to get filled. So the way of calculating
% fast rate az and el without pointing model which was used in QUaD
% doesn't work.

% tracker.horiz_mount does not include encoder zero points but
% pmac.fast_az_pos etc does - add them to the former
horiz_mount=d.antenna0.tracker.horiz_mount*180/pi+double(d.antenna0.tracker.encoder_off)/3.6e6;

% horiz_mount does include scan and user offsets - remove them
horiz_mount=horiz_mount-double(d.antenna0.tracker.scan_off)/3.6e6;
horiz_mount=horiz_mount-d.antenna0.tracker.horiz_off*180/pi;

% expand to fast rate
samprate=length(d.antenna0.pmac.fast_az_pos)/size(d.antenna0.tracker.horiz_mount,1);
horiz_mount_az=repmat(horiz_mount(:,1)',samprate,1);
horiz_mount_az=cvec(horiz_mount_az);
horiz_mount_el=repmat(horiz_mount(:,2)',samprate,1);
horiz_mount_el=cvec(horiz_mount_el);

% get encoder conversion factors (DK changed from bicep1 to bicep2)
enc_mul= double(d.antenna0.tracker.encoder_mul(1,:))./360;

% find actual x/y offset at each time step
d.eloff=double(d.antenna0.pmac.fast_el_pos)/enc_mul(2)-horiz_mount_el;
d.azoff=double(d.antenna0.pmac.fast_az_pos)/enc_mul(1)-horiz_mount_az;
d.azoff=mod(d.azoff+180,360)-180;

% rawpoint is just use crude conversion from fast rate register
d.rawpoint.hor.az=double(d.antenna0.pmac.fast_az_pos)/enc_mul(1)-cvec(repmat((double(d.antenna0.tracker.encoder_off(:,1))/3.6e6)',samprate,1));
d.rawpoint.hor.el=double(d.antenna0.pmac.fast_el_pos)/enc_mul(2)-cvec(repmat((double(d.antenna0.tracker.encoder_off(:,2))/3.6e6)',samprate,1));
d.rawpoint.hor.dk=double(d.antenna0.pmac.fast_dk_pos)/enc_mul(3)-cvec(repmat((double(d.antenna0.tracker.encoder_off(:,3))/3.6e6)',samprate,1));

return
