function supfac=calc_supfac(s)
% supfac=calc_supfac(s)
%
% Calculate filter/beam supresion factor by comparing the expectation
% value of a set of signal only sims to the mean value obtained in
% each aps bin

 
% calc filter/beam supression factor
for j=1:length(s)
  % calc supression factors by comparing signal only sim output to the
  % expectation value in each bin
  supfac(j).tt=s(j).mean(:,1)./s(j).expv(:,1);
  supfac(j).ee=s(j).mean(:,3)./s(j).expv(:,3);
  supfac(j).bb=s(j).mean(:,4)./s(j).expv(:,4);
  supfac(j).l=s(j).l;
  
  % calc error on supfac due to error on mean of signal only
  % simulations (due to finite number of sims)
  supfac(j).tteu=(s(j).mean(:,1)+s(j).eom(:,1))./s(j).expv(:,1)-supfac(j).tt;
  supfac(j).eeeu=(s(j).mean(:,3)+s(j).eom(:,3))./s(j).expv(:,3)-supfac(j).ee;
  supfac(j).bbeu=(s(j).mean(:,4)+s(j).eom(:,4))./s(j).expv(:,4)-supfac(j).bb;
  supfac(j).ttel=supfac(j).tt-(s(j).mean(:,1)-s(j).eom(:,1))./s(j).expv(:,1);
  supfac(j).eeel=supfac(j).ee-(s(j).mean(:,3)-s(j).eom(:,3))./s(j).expv(:,3);
  supfac(j).bbel=supfac(j).bb-(s(j).mean(:,4)-s(j).eom(:,4))./s(j).expv(:,4);
  
  % calc supfac from TE - this is messy due to zero crossings
  supfac(j).te=s(j).mean(:,2)./s(j).expv(:,2);
  supfac(j).teeu=(s(j).mean(:,2)+s(j).eom(:,2))./s(j).expv(:,2)-supfac(j).te;
  supfac(j).teel=supfac(j).te-(s(j).mean(:,2)-s(j).eom(:,2))./s(j).expv(:,2);
 
  % for freq/rx cross calc ET sup fac as well
  if(size(s(j).Cs_l,2)==9)
    supfac(j).et=s(j).mean(:,7)./s(j).expv(:,7);
    supfac(j).eteu=(s(j).mean(:,7)+s(j).eom(:,7))./s(j).expv(:,7)-supfac(j).et;
    supfac(j).etel=supfac(j).et-(s(j).mean(:,7)-s(j).eom(:,7))./s(j).expv(:,7);
  end

end

% attempt to fit supression factor to smooth it out
if(0)
  for j=1:length(s)
    x=supfac(j).l;
    y=supfac(j).tt; e=mean([supfac(j).tteu,supfac(j).ttel],2);
    %supfac(j).ptt=matmin('chisq',[1,zeros(1,6)],[],'polyval2',y,[],x);
    supfac(j).ptt=polyfit(x,y,5);
    
    y=supfac(j).ee; e=mean([supfac(j).eeeu,supfac(j).eeel],2);
    %supfac(j).pee=matmin('chisq',[1,zeros(1,6)],[],'polyval2',y,[],x);
    supfac(j).pee=polyfit(x,y,5);
  end
end
  
if(0)
  % replace supfac with fit val
  for j=1:length(s)
    supfac(j).tt=polyval(supfac(j).ptt,supfac(j).l);
    supfac(j).ee=polyval(supfac(j).pee,supfac(j).l);
  end
end
  
for j=1:length(s)
  % to avoid TE zero crossing problem
  % assume cross spectra supressed by geometric mean of corresponding
  % auto spectra
  % - this approx does not appear to work that well...
  % In as much as supfac really differ between spectra at a given freq then
  % the above is not quite right.
  % For T100E150 supfac should be gmean(TT100,EE150) whereas for
  % E100T150 the supfac should be gmean(EE100,TT150)
  supfac(j).te=sqrt(supfac(j).tt.*supfac(j).ee);
  supfac(j).tb=sqrt(supfac(j).tt.*supfac(j).bb);
  supfac(j).eb=sqrt(supfac(j).ee.*supfac(j).bb);
end

% construct the multiplier to be applied to obs spectra
for i=1:length(s)
  % the multiplier is the reciprocal of the window function - hence rwf
  tt=supfac(i).tt; ee=supfac(i).ee; bb=supfac(i).bb;
  te=supfac(i).te; tb=supfac(i).tb; eb=supfac(i).eb;
  if(size(s(i).Cs_l,2)==6)
    supfac(i).rwf=1./[tt,te,ee,bb,tb,eb];
  else
    supfac(i).rwf=1./[tt,te,ee,bb,tb,eb,te,tb,eb];
  end
end

return
