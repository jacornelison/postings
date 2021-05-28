function scratch()

%% addpaths
addpath('Z:/pipeline/util/')
addpath('Z:/pipeline/external_lib/')

%% Load the things

load('../data/b3rpsfiles.mat')
f = fopen('../data/20180116_RPS_monitor_3.csv','r');

data = textscan(f,'%f%f%f%f%f%f%f%f%f%s','delimiter',',','HeaderLines',1);
%%
% convert python time.time() into mjd
% time.time() epoch is 01 Jan 1970 UTC
time = data{1}/3600/24 + date2mjd(1970,01,01,0,0,0);

%%
%time, srctemp, PID, tiltout+, tiltout-,tilttemp+,tilttemp-
az_offs = [-4.5,-13,4,-4.5,4.0];
dk = 0;
for schind = 1:length(sch)
    
    if dk ~= sch{schind}.dk_ofs
        count = 1;
        dk = sch{schind}.dk_ofs;
    else
        count = count+1;
    end
    sch{schind}.az_offs = az_offs(count);
end


pm = [];
pm.meas_name = {'170203 02:32:23'};
pm.mjd = 5.7787e+04;
pm.flex_cos = 0;
pm.flex_sin = 0;
pm.az_tilt_ha = 0.0046;
pm.az_tilt_lat = 0.0044;
pm.el_tilt = 0.0286;
pm.collim_x = 0;
pm.collim_y = 0;
pm.collim_mag = 0;%-0.5802;
pm.collim_dir = 0;%143.9816;
pm.az_zero = 0;
pm.el_zero = 0;
pm.err_flex_cos = 0;
pm.err_flex_sin = 0;
pm.err_az_tilt_ha = 0.0614;
pm.err_az_tilt_lat = 0.0610;
pm.err_el_tilt = 0.2295;
pm.err_collim_x = 0;
pm.err_collim_y = 0;
pm.err_collim_mag = 0.0585;
pm.err_collim_dir = 5.8488;
pm.err_az_zero = 0.4323;
pm.err_el_zero = 0.0604;
pm.rms_res = 0.0050;
pm.numpts = 32;



%% plot a scanset
%close all
figure(1)
clf; hold on;

figure(2)
clf; hold on;


figure(3)
clf; hold on;

map = colormap;
lim = [-1,1]*0.2;
T0 = 0.105*1000;
Tcal = 0.04;
count = 0;
Tthresh = 0.2;
Vthresh = 0.01;
DtiltT = (data{6}-data{7})*100;
DtiltV = (data{4}-data{5})*1000/1666;
DT = DtiltT;
DV = DtiltV;
for schind = 1:3%1:length(sch)

    nrows = sch{schind}.nrows;
    nrps = sch{schind}.nrps;
    [Toff,Voff,Tm,Vm] = deal(zeros(1,nrows*nrps));
    
    
    for rowind = 1:nrows
        for rpsind = 1:nrps
            scanind = (rowind-1)*nrps+rpsind;
            scan = sch{schind}.scans(scanind);
            t1 = scan.t1;
            t2 = scan.t2;
            tind = find(time > t1 & time < t2);
            if ~isempty(tind)
                
                Ttilt = DT(tind);
                Vtilt = DV(tind);
                Tm(scanind) = nanmean(Ttilt);
                Vm(scanind) = nanmean(Vtilt);
                
                
                if scanind==1
                    tstart = time(tind(1));
                    
                else
                    dT = diff(Tm([-1,0]+scanind));
                    dV = diff(Vm([-1,0]+scanind));
                    if abs(dT)>Tthresh | abs(dV)>Vthresh
                        DT = DT-dT;
                        DV = DV-dV;
                        Ttilt = DT(tind);
                        Vtilt = DV(tind);
                        Tm(scanind) = nanmean(Ttilt);
                        Vm(scanind) = nanmean(Vtilt);
                    end
                end
                
                
                
                tscan = (time(tind)-tstart)*24*60;
                                
                
%                 figure(1)
%                 subplot(2,1,1)
%                 hold on
%                 plot(tscan,DtiltV(tind))
%                 grid on
%                 subplot(2,1,2)
%                 hold on
%                 plot(tscan,Vtilt)
%                 grid on
%                 
%                 
%                 figure(2)
%                 subplot(2,1,1)
%                 hold on
%                 plot(tscan,DtiltT(tind))
%                 grid on
%                 subplot(2,1,2)
%                 hold on
%                 plot(tscan,Ttilt)
%                 grid on
%                 
            end
        end
                az = [0,9,9,0]+sch{schind}.az_offs+2.5;
                el = 90-[-1 -1 1 1]+scan.el_ofs-2;
                dk = [1 1 1 1]*scan.dk_ofs-1.5;
                [x,y,phi] = beam_map_pointing_model(az,el,dk,pm,'bicep3',rpsopt.mirror,rpsopt.source,[]);
                figure(3)
                mV = mean(Vm((rowind-1)*nrps+[1:13]));
                fill(x,y,mV)
                xlim([-25,25])
                ylim([-1 1]*25)
                            
        
    end
    
end



figure(3)
prx = 2 * sind(p.r / 2) .* cosd(p.theta) * 180 / pi;
pry = 2 * sind(p.r / 2) .* sind(p.theta) * 180 / pi;
plot(prx,pry,'rx')
grid on
colorbar()

%%

figure(3)
plot((data{4}-data{5})*1000/1666,(data{6}-data{7}),'.')


