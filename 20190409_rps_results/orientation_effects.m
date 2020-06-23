function orientation_effects(a,b,offs,pol,C)
% function orientation_effects(a,b,offs,pol,C)


%% Investigate the systematic bias on psi due to RPS orientation

% Our model from Colin's posting
% A(theta) = [(cosd(psi)*cosd(b)+sind(psi)*sind(a)*sind(b))*cosd(th)+sind(psi)*cosd(a)*sind(th)].^2

% Initialize parameters
th = -180:30:180;
sig = 0.3;
sig2 = 0.1;
%offs = -2;
dk = [0 45 90 135];
freepar.free = [1 1 1 1 1];
freepar.lb = [-360 -1 -1e4 -1e4 -1e4];
freepar.ub = [360 1 1e4 1e4 1e4];

siga = [0.17 0.33 0.08 0.24 0.19 0.14 0.13 0.26 0.24 0.16 0.24 0.29 0.13 0.24 0.11 0.08 0.15 0.20 0.05 0.17];
%sigb = [0.17 0.33 0.08 0.24 0.19 0.14 0.13 0.26 0.24 0.16 0.24 0.29 0.13 0.24 0.11 0.08 0.15 0.20 0.05 0.17];

%sigb = [];
% pol = 'a'
% a = -10;
% b = 10;
parm_all = [];
parm_allx = [];

tiles = normrnd(offs,sig,20,1);
chans = 64;
for j = 1:4
    dks = dk(j);
    parm_psi = [];
    parm_x = [];
    for k = 1:length(tiles)
    for i = 1:chans
        
        psi = dks+tiles(k)+normrnd(0,siga(k));
        A = ((cosd(psi)*cosd(b)+sind(psi)*sind(a)*sind(b))*cosd(th)+sind(psi)*cosd(a)*sind(th)).^2;
        guess = [psi 0.000 0 0 1/2];
        
        % Estimate parameters
        [aparam, aerr, agof, astat, acov] = matmin('chisq',...
            guess, freepar,	'rps_get_mod_model',A,1,th);
        
        parm_psi = [parm_psi aparam(1)-dks];
        parm_x = [parm_x aparam(2)];
    end
    end
    parm_all = [parm_all; parm_psi];
    parm_allx = [parm_allx; parm_x];
end

a1 = parm_all(1,:);
a2 = parm_all(2,:);
a3 = parm_all(3,:);
a4 = parm_all(4,:);

B = 1;
C = C+offs;
clr = [0,0,0.7];
figfold = 'figs3/';

fig = figure(2)
set(fig,'Position',[500 0 [200 300]*1.7])
pname = ['postings/20180731_rps_something/' figfold sprintf(['rps_orientation_a_%02i_b_%02i_offs_%02i_pol_' pol],a,b,offs)];
clf(fig)

subplot(3,2,1)
plot(a1,a2,mk,'Color',clr)
hold on
c = polyfit(a1,a2,1);
plot(([-10 10]+C),([-10 10]+C)*c(1)+c(2),'--','Color',clr)
plot(([-10 10]+C),([-10 10]+C),'k--','LineWidth',1)
title(sprintf('DK %2.0f VS DK %2.0f',dk(1),dk(2)))
xlabel('Phi (^o)')
ylabel('Phi (^o)')
ylim([-4 4]*B+C)
xlim([-4 4]*B+C)
axis square
grid on

subplot(3,2,2)
plot(a1,a3,mk,'Color',clr)
hold on
c = polyfit(a1,a3,1);
plot(([-10 10]+C),([-10 10]+C)*c(1)+c(2),'--','Color',clr)
plot(([-10 10]+C),([-10 10]+C),'k--','LineWidth',1)
title(sprintf('DK %2.0f VS DK %2.0f',dk(1),dk(3)))
xlabel('Phi (^o)')
ylabel('Phi (^o)')
ylim([-4 4]*B+C)
xlim([-4 4]*B+C)
axis square
grid on

subplot(3,2,3)
plot(a1,a4,mk,'Color',clr)
hold on
c = polyfit(a1,a4,1);
plot(([-10 10]+C),([-10 10]+C)*c(1)+c(2),'--','Color',clr)
plot(([-10 10]+C),([-10 10]+C),'k--','LineWidth',1)
title(sprintf('DK %2.0f VS DK %2.0f',dk(1),dk(4)))
xlabel('Phi (^o)')
ylabel('Phi (^o)')
ylim([-4 4]*B+C)
xlim([-4 4]*B+C)
axis square
grid on

subplot(3,2,4)
plot(a2,a3,mk,'Color',clr)
hold on
c = polyfit(a2,a3,1);
plot(([-10 10]+C),([-10 10]+C)*c(1)+c(2),'--','Color',clr)
plot(([-10 10]+C),([-10 10]+C),'k--','LineWidth',1)
title(sprintf('DK %2.0f VS DK %2.0f',dk(2),dk(3)))
xlabel('Phi (^o)')
ylabel('Phi (^o)')
ylim([-4 4]*B+C)
xlim([-4 4]*B+C)
axis square
grid on

subplot(3,2,5)
plot(a2,a4,mk,'Color',clr)
hold on
c = polyfit(a2,a4,1);
plot(([-10 10]+C),([-10 10]+C)*c(1)+c(2),'--','Color',clr)
plot(([-10 10]+C),([-10 10]+C),'k--','LineWidth',1)
title(sprintf('DK %2.0f VS DK %2.0f',dk(2),dk(4)))
xlabel('Phi (^o)')
ylabel('Phi (^o)')
ylim([-4 4]*B+C)
xlim([-4 4]*B+C)
axis square
grid on

subplot(3,2,6)
plot(a3,a4,mk,'Color',clr)
hold on
c = polyfit(a3,a4,1);
plot(([-10 10]+C),([-10 10]+C)*c(1)+c(2),'--','Color',clr)
plot(([-10 10]+C),([-10 10]+C),'k--','LineWidth',1)
xlabel('Phi (^o)')
ylabel('Phi (^o)')
title(sprintf('DK %2.0f VS DK %2.0f',dk(3),dk(4)))
ylim([-4 4]*B+C)
xlim([-4 4]*B+C)
axis square
grid on

saveas(fig,pname,'png')
