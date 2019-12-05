%% Attenuation Curves
load('attenuation_curves_03dec2019.mat')

for i = 1:2
    fig = figure();
    if i == 1
        x = hughes(:,1);
        A = hughes(:,2);
        lbl = 'Hughes';
    else
        x = milli(:,1);
        A = milli(:,2);
        lbl = 'Millitech';
    end
    A = 10*log10(A/max(A));
    plot(x,A)
    grid on
    ylim([-40,5])
    title([lbl ' Attenuation Curve'])
    xlabel('Dial Setting')
    ylabel('Amplitude (dB)')
end

%% Spectra
load('converted_spectra.mat')

[f1, ind] = sort(new_spect_2040(:,1));
A1 = new_spect_2040(ind,2);
[f2, ind] = sort(new_spect(:,1));
A2 = new_spect(ind,2);

A = [A1(f1<24);A2];
f = [f1(f1<24);f2];

fig = figure();
clf(fig)
plot(0,0,'b',0,0,'b--',0,0,'r',0,0,'r--')
hold on
plot(f,A,'k')

bc = [30.8, 41.5];
bw = [0.235 0.28];

plot([1 1]*bc(1),[-50,0],'b')
plot([1 1]*bc(1)*(1-bw(1)/2),[-50,0],'b--')
plot([1 1]*bc(1)*(1+bw(1)/2),[-50,0],'b--')

plot([1 1]*bc(2),[-50,0],'r')
plot([1 1]*bc(2)*(1-bw(2)/2),[-50,0],'r--')
plot([1 1]*bc(2)*(1+bw(2)/2),[-50,0],'r--')

ylim([-42,-25])
xlim([20, 50])
title('30/40GHz BBNS Spectrum')
xlabel('Frequency (GHz)')
ylabel('Amplitude (dBm)')
legend('30GHz BC','30GHz BW','40GHz BC','40GHz BW')
grid on

%% Temperature Dependence
load('tempvsvoltage_table_02dec2019.mat')

ind = 5000; % don't include warm-up period
T = TVtable{ind:end,2}*1000;
V = TVtable{ind:end,3};
V = V/median(V);

% Get the thermal dependence
c = polyfit(T,V,1);
M = c(1)*T+c(2);

fig = figure();
subplot(1,1,1)
p = plot(T,V,T,M);
p(2).LineWidth = 2;
xlabel('Temperature (K)')
ylabel('Amplitude (Normalized)')
text(310,1.02,sprintf('Slope: %0.2f %%/K',c(1)*100))
legend('Data','Model')
grid on

%% Beams
load('beam_maps_03dec2019.mat')

fig = figure();
clf(fig)
fig.Position = [200 300 800 400];
for i = 1:2
    if i==1
        map = emap;
        mlab = "E-plane";
        k = -1;
    else
        map = hmap;
        mlab = "H-plane";
        k = 0;
    end
    
    theta = map(:,1);
    amp = map(:,2);
    amp = amp/max(amp);
    
    fun = @(p) nansum((amp-ampgaussmf(p,theta)).^2);
    param = fminsearch(fun,[1 18 0 0]);
    
    theta = theta-param(3);
    param(3) = 0;
    
    subplot(2,2,1+k+i)
    plot(theta,ampgaussmf(param,theta))
    hold on
    plot(theta,amp,'k.')
    grid on
    title([mlab " Data vs Model"])
    legend("Model","Data")
    xlabel("Theta (deg)")
    ylabel("Amplitude (uV)")
    ylim([0,1.1])
    
    text(-55,0.75,['$\sigma$' sprintf(' = %0.1f',param(2)) '$^\circ$'],'Interpreter','latex')
    
    
    subplot(2,2,2+k+i)
    plot(theta,(amp-ampgaussmf(param,theta))/max(amp))
    grid on
    title('Best-Fit Fractional Residuals')
    xlabel("Theta (deg)")
    ylim([-1,1]*0.02)
        
end


% Subfunction for beam fitting
function amp = ampgaussmf(p,x)
if length(p)>3
    amp = p(1)*gaussmf(x,[p(2) p(3)])+p(4);
else
    amp = p(1)*gaussmf(x,[p(2) p(3)]);
end

end