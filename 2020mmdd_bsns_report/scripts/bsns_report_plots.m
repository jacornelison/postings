function bsns_report_plots()

%% Global Initializations
close all
pltind = 1;

%% Total power
bsns_label = {'30 GHz','90 GHz','150 GHz','220 GHz'};
bsns_power = [12.6, 2.3, -0.7, -40]; % dBm
src_color = {[0, 0.4470, 0.7410], [0.6350, 0.0780, 0.1840], [0.4660, 0.6740, 0.1880],	[0.4940, 0.1840, 0.5560]};
%% Dynamic Range
% Attenuation curves in format [dial setting, amplitude (dB); ...]

att_names = {...
    'milli','hughes';...
    'rocky','severe_backlash';...
    'att_2','att_3';...
    'att_1','att_2'};

att_titles = {...
    '30 GHz Att 1','30 GHz Att 2';...
    '90 GHz Att 1','90 GHz Att 2';...
    '150 GHz Att 1','150 GHz Att 2';...
    '220 GHz Att 1','220 GHz Att 2'};


figure(pltind); pltind=pltind+1;
clf
set(gcf,'Position',[50,50,700,700])
for srcind = 1:4
    att_curves = load(['data\' bsns_label{srcind} '\attenuation_curves.mat']);
    

    
    for attind = 1:2
        att_obj = att_curves.(att_names{srcind,attind});
        
        %sortind = 1:length(att_obj(:,1));
        [B, sortind] = sort(att_obj(:,1));
        
            % For 90 GHz, downsample.
        if srcind == 2 | srcind == 3
            att_obj = att_obj(sortind,:);
            [uval,uind] = unique(att_obj(:,1));
            m = max(att_obj(:,1));
            dial_int = 0:m/10:m;
            pwr_int = interp1(att_obj(uind,1),att_obj(uind,2),dial_int);
            att_obj = dial_int';
            att_obj(:,2) = pwr_int';
            [B, sortind] = sort(att_obj(:,1));
        end
        
        subplot(4,2,2*srcind+attind-2)
        plot(att_obj(sortind,1),att_obj(sortind,2),'-o','Color',src_color{srcind},...
        'MarkerFaceColor',src_color{srcind});
        grid on
        title(att_titles{srcind,attind})
        %ylim([nanmin(att_obj(:,2))*1.1,5])
        ylim([-55,5])
        
        
        yticks(-45:10:5)
        if srcind == 1
            set(gca,'XDir','reverse')
        elseif srcind == 4
            xlabel('Dial Setting')
        end
        
        if attind ==1
            ylabel('Amplitude (dB)')
        end
    end
    
end


%% Output Stability - Room Temp

figure(pltind); pltind=pltind+1;
clf
set(gcf,'Position',[300,200,900,400])
pltlegend = {'30 GHz (South Pole, 250K Ambient)','90 GHz (In-Lab,300K Ambient)','','220 GHz (In-Lab,300K Ambient)'};
src_order = [2 4 1];

for srcind = src_order
    for pind = 1:2
    
    subplot(1,2,pind)
    src_adev = load(['data\' bsns_label{srcind} '\bsns_src_stability_final.mat']);
    set(gca, 'YScale', 'log','XScale','log')
    if pind == 2
            
        errorbar(src_adev.Ptau,(src_adev.Pad),(src_adev.Pade),'-o','Color',src_color{srcind},...
        'MarkerFaceColor',src_color{srcind},'CapSize',10)
        ylabel('\sigma_P(\tau) Fractional Units')
        title('Output Stability')
    else 
        errorbar(src_adev.Ttau,(src_adev.Tad),(src_adev.Tade),'-o','Color',src_color{srcind},...
        'MarkerFaceColor',src_color{srcind},'CapSize',10)
        ylabel('\sigma_T(\tau) Fractional Units')
        title('Temperature Stability')
    end
    hold on
    
    xlabel('\tau (seconds)')
    ylim([5e-10,7e-1])
    xlim([1e-1,5e4])
    grid on  
    end
end

pltlegend = pltlegend(src_order);
pltlegend{end+1} = 'Target Stability';
subplot(1,2,2)
plot([1e-1, 4e4],1e-2*[1, 1],'k--')
legend(pltlegend)

subplot(1,2,1)
plot([1e-1, 4e4],exp(-0.5*[1e-1, 4e4]/1e-1),'k-')

%% Output Stability - Timestreams
if 1
figure(pltind); pltind=pltind+1;
clf
set(gcf,'Position',[300,200,900,400])
pltlegend = {'30 GHz (South Pole, 250K Ambient)','90 GHz (In-Lab,300K Ambient)','','220 GHz (In-Lab,300K Ambient)'};
src_order = [2 4 1];

for srcind = src_order
    src_adev = load(['data\' bsns_label{srcind} '\bsns_src_stability_final.mat']);
        
    for pind = 1:2
    
    subplot(1,2,pind)
    %set(gca, 'YScale', 'log','XScale','log')
    if pind == 2
            
        plot(src_adev.tint,src_adev.Pint+1,'Color',src_color{srcind})
        ylabel('Median-Normalized Amplitude')
        title('Output Timestream')
        ylim([0.92, 1.12])
    else 
        plot(src_adev.tint,src_adev.Tint/median(src_adev.Tint),'Color',src_color{srcind})
        ylabel('Median-Normalized Temperature')
        title('Temperature Timestream')
        ylim([0.985, 1.02])
    end
    hold on
    
    xlabel('Time (seconds)')
    %ylim([5e-4,2])
    %xlim([1e-1,5e4])
    grid on  
    end
end

legend(pltlegend(src_order))
end
%% Bandpass
% Spectra are in 'data/# GHz/final_spectra.mat'
% Format f = frequency in GHz and A is amplitude in dB (peak-normalized)

% 3040 GHz - from spectrum analyzer
load('data\30 GHz\final_spectra.mat')

% 90 GHz

% 150 GHz

% 220 GHz

%%%%%%%%%%%%%%%%%%%
%% RPS
%%%%%%%%%%%%%%%%%%%
%% Pol efficiency / Repeatability
% Need to do per-freq since we use the same wire grid for each.

% 3040 GHz

% 90 GHz

% 150 GHz

% 220 GHz

%% Collimation Stability


%% Gain VS. Rotation Angle
% Need to do per-freq since we use different chassis / wiring

% 3040 GHz

% 90 GHz

% 150 GHz

% 220 GHz

%% Beams
% This is beams through the wire grid.

% 3040 GHz

% 90 GHz

% 150 GHz

% 220 GHz



