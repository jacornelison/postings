function [rp, thp, psii] = pointing_checks(az,el,dk,p,ch,rpsopt)

mount = rpsopt.mount;
mirror = rpsopt.mirror;
source = rpsopt.source;
source.distance = 1e7;
source.height = source.distance*tan(source.el);

[rp,thp,psii] = deal(NaN(size(ch)));

for i=1:length(ch)
    if isfield(p, 'expt')
        p = rmfield(p, 'expt');
    end
    pind = structcut(p, ch(i));
    %[rp(i),thp(i),psii(i)] = keck_beam_map_pointing(az(i),el(i),dk(i),mount,mirror,source,pind,'NoMirror',NoSource);
    [rp(i),thp(i),psii(i)] = keck_beam_map_pointing(source.azimuth,source.el,dk(i),mount,mirror,source,pind,'NoMirror','NoSource');
end

xp = 2 * sind(rp / 2) .* cosd(thp) * 180 / pi;
yp = 2 * sind(rp / 2) .* sind(thp) * 180 / pi;
