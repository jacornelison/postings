function fd = rps_cut_fit_data(fd,p)
% function fd = rps_cut_fit_data(fd)
% Cut parameters are currently hard coded.
% go back and fix this.

fprintf('Cutting wonky channels\n.\n.\n.\n')



prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;

dk = unique(fd.dk);
cut_sch = [08 06 07 07 07 07 08 09 09 05 06 06 04 04 05 05 05 09 11];
cut_row = [11 14 18 02 11 17 17 04 11 11 13 04 17 09 07 10 12 15 02];

cut_raster = zeros(size(fd.ch));
for i = 1:length(cut_sch)
    cut_raster = cut_raster | (fd.sch==cut_sch(i) & fd.row==cut_row(i));
end

cut_ind = ...
    abs(prx(fd.ch)-fd.bparam(:,1))>1 |...
    abs(pry(fd.ch)-fd.bparam(:,2))>1 |...
    (fd.phi_d < -4 & fd.phi_d > 02 & fd.inda == 1)|... % No pol a's should be 90 degrees.
    (fd.phi_d < 86 & fd.phi_d > 92 & fd.inda ~= 1)|... % No pol b's should be 0 degrees.
    abs(fd.aparam(:,2)) > 0.05 |...
    any(abs(fd.res_perc)>0.1,2)|...
    abs(prx(fd.ch)-fd.bparam(:,1))>0.2 |...
    abs(pry(fd.ch)-fd.bparam(:,2))>0.2 |...
    abs(fd.aparam(:,3))>0.03 |... % High miscollimation indicates poor source stability
    abs(fd.aparam(:,4))>0.03 |...
    log10(fd.net)>3; % Cut noisy detectors

fd = structcut(fd,~cut_ind & ~cut_raster);


cutchans_a = unique(fd.ch(fd.inda==1));
cutchans_b = unique(fd.ch(fd.inda~=1));
if 1
    % Use only data with pairs
    inda = zeros(size(cutchans_a));
    indb = zeros(size(cutchans_b));
    for k = 1:20
        for l = 1:8
            for mi = 1:8
                cha = find(p.tile==k & p.det_row==l & p.det_col==mi & strcmp(p.pol,'A'));
                chb = find(p.tile==k & p.det_row==l & p.det_col==mi & strcmp(p.pol,'B'));
                
                if (~isempty(cha) & ~isempty(chb)) & p.mce(cha)==0
                    if  any(cutchans_a==chb) & any(cutchans_b==cha)
                        inda(cutchans_a==chb) = 1;
                        indb(cutchans_b==cha) = 1;
                    end
                elseif (~isempty(cha) & ~isempty(chb)) & p.mce(cha)~=0
                    if  any(cutchans_a==cha) & any(cutchans_b==chb)
                        inda(cutchans_a==cha) = 1;
                        indb(cutchans_b==chb) = 1;
                    end
                end
            end
        end
    end
    
    use_ind = ismember(fd.ch,cutchans_a(inda==1)) | ismember(fd.ch,cutchans_b(indb==1));
    fd = structcut(fd,use_ind);
end
