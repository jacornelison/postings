function plot_bmspec(aps,leg,spectype,lmax,freq,jack,dp,rho)
% function plot_bmspec(aps,leg,spectype,lmax,freq,jack,dp,rho)
% 
% Utility function for plotting beam map spectra
% Does one of the 6-panel plots at a time
% Can send in as many aps as you like - just loops over length of aps
% 
% leg is a char array with the names of the spectra for the legend 
% 
% lmax/spectype control the axes of the plot
% freq/jack/dp just control the title

if ~isempty(rho)
  if length(rho) ~= length(aps)
    disp('Rho not same length as aps, not plotting')
    plotrho = 0;
  else
    plotrho = 1;
  end
else
  plotrho = 0;
end

switch spectype
  case 'TT'
    specswitch = 1;
  case 'TE'
    specswitch = 2;
  case 'EE'
    specswitch = 3;
  case 'BB'
    specswitch = 4;
    % Model curves - r = 0.1
    specname = 'input_maps/official_cl/camb_planck2013_r0p1_noE.fits';
    r0p1 = load_cmbfast(specname);
    % Lensing
    specname = 'input_maps/official_cl/camb_planck2013_r0_lensing.fits';
    lens = load_cmbfast(specname);
  case 'TB'
    specswitch = 5;
  case 'EB'
    specswitch = 6;
end

switch lmax
  case '200'
    xra = [0 200];
  case '500'
    xra = [0 500];
end

% Change this if we want different colors, etc...
linestyles{1} = '-b';
linestyles{2} = '-g';
linestyles{3} = '-m';

% Make the plot
switch spectype
  case 'BB'
    % Add BB - r = 0.01 and r = 0.2 + lensing
    semilogy(r0p1.l,0.1*r0p1.Cs_l(:,4),'--r');
    hold on;
    semilogy(r0p1.l,2*r0p1.Cs_l(:,4) + lens.Cs_l(:,4),'--r');
    for ii = 1:length(aps)
      semilogy(aps(ii).l(2:end),...
               aps(ii).Cs_l(2:end,specswitch),...
               linestyles{ii},'LineWidth',2);
      hold on;
      if plotrho
        text(0.65*xra(2),2e-5,...
             ['\rho = ' sprintf('%0.4f',rho(ii))]);
      end
    end
    
    xlabel('ell');
    ylabel('l(l+1)C_l/2\pi (uK^2)');

  otherwise
    for ii = 1:length(aps)
      plot(aps(ii).l(2:end),...
           aps(ii).Cs_l(2:end,specswitch),...
           linestyles{ii},'LineWidth',2); 
      hold on;
      switch spectype
        case 'EE'
          legend(leg,'Location','NorthWest')
      end
    end
end

% Tailor limits to type of spectrum

xlim(xra);

switch spectype
  case 'TT'
    ylim([0 6500])
  case 'TE'
    ylim([-50 50])
    set(gca,'YTick',[-50 -25 0 25 50]);
  case 'EE'
  case 'BB'
    ylim([1e-5 1e0])
    set(gca,'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0]);
  case 'TB'
    ylim([-10 10]);
    set(gca,'YTick',[-10 -5 0 5 10]);
  case 'EB'
    ylim([-0.5 0.5]);
    set(gca,'YTick',[-0.5 -0.25 0 0.25 0.5]);
end

grid on;
set(gca,'YMinorGrid','off')
title(sprintf('%s %s: jack%s dp%s',freq,spectype,jack,dp))


return

