m1 = load('out/B3_SP2018_20181217_run10_OE/inter_results.mat');
m2 = load('out/B3_SP2019_20191228_run11_OE/inter_results.mat');
mk = {[0,0,0.7],[0.7,0,0],[0,0.7,0.2],[0.5,0,0.5]};
p = get_array_info(20191228);
figure(1)
fpos = [1000,300,800,800];
paspect = 5;
set(gcf,'Position',fpos);
set(gcf,'PaperPosition',[0,0,paspect,paspect*fpos(4)/fpos(3)]);
kboltz = 1.38e-23;
dnu = 25e9;

for k =1:3
	if k==1
		dp1 = m1.PSat{1}*1e12;
		dp2 = m2.PSat{1}*1e12;
		figfile = '~/20191230_bicep3_opteff/figs/psat_compare_run10_to_run11';
		lim = [0,350];
		lim2 = [0,5];
		stp = 70;
		t1 = 'Run 11 P_{sat} (pW)';
		t2 = 'Run 10 P_{sat} (pW)';
		t0 = 'Run 11 / Run 10 P_{Sat} Comparisons';
		t3 = 'P_{Sat,10}-P_{Sat,11} (pW)';
	elseif k==2
		dp1 = m1.dPdT*1e12;
		dp2 = m2.dPdT*1e12;
		figfile = '~/20191230_bicep3_opteff/figs/dpdt_compare_run10_to_run11';
		lim = [0,0.25];
		lim2 = [0,0.05];
		stp = 70;
		t1 = 'Run 11 dP/dT (pW/K)';
		t2 = 'Run 10 dP/dT (pW/K)';
		t0 = 'Run 11 / Run 10 dP/dT Comparisons';
		t3 = 'dPdT_{10}-dPdT_{11} (pW)';
	else
		dp1 = m1.dPdT/kboltz/dnu;
		dp2 = m2.dPdT/kboltz/dnu;
		figfile = '~/20191230_bicep3_opteff/figs/oe_compare_run10_to_run11';
		lim = [0,0.30]+0.15;
		lim2 = [0,0.1];
		stp = 70;
		t1 = 'Run 11 OE';
		t2 = 'Run 10 OE';
		t0 = 'Run 11 / Run 10 OE Comparisons';
		t3 = 'OE_{10}-OE_{11}';
	end

	for i = 1:4
	subplot(2,2,i)
	clf
	end

	edg = lim(1):diff(lim)/stp:lim(2);
	edg2 = lim2(1):diff(lim2)/stp:lim2(2);
	for i=1:length(unique(p.mce))
		ind = p.mce==(i-1);
		subplot(2,2,1)
		hold on
		N = histc(dp2(ind),edg);
		bar(edg,N,'FaceColor',mk{i})

		subplot(2,2,2)
		hold on
		plot(dp1(ind),dp2(ind),'.','Color',mk{i})

		subplot(2,2,3)
		hold on
		pdiff = dp1-dp2;
		pdiff = pdiff(p.mce==(i-1)&inrange(dp1,lim(1),lim(2))&inrange(dp1,lim(1),lim(2)));
		bar(edg2,histc(pdiff,edg2),'FaceColor',mk{i})
		
		subplot(2,2,4)
		hold on
		N = histc(dp1(ind),edg);
		bar(edg,N,'FaceColor',mk{i})
	end

	subplot(2,2,1)
	title(t0)
	xlim(lim)
	ylabel('N')
	xlabel(t1)
	grid on
	legend('MCE0','MCE1','MCE2','MCE3')

	subplot(2,2,2)
	plot([-1 1]*1000,[-1 1]*1000,'k--')
	xlim(lim)
	ylim(lim)
	legend('MCE0','MCE1','MCE2','MCE3')
	title('')
	grid on
	xlabel(t2)
	ylabel(t1)

	subplot(2,2,3)
	xlim(lim2)
	ylabel('N')
	xlabel(t3)
	grid on
	legend('MCE0','MCE1','MCE2','MCE3')

	subplot(2,2,4)
	xlim(lim)
	ylabel('N')
	xlabel(t2)
	grid on
	legend('MCE0','MCE1','MCE2','MCE3')

	print('-dpng',figfile)

end


