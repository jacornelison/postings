function [noise mn]=estimatenoise(mapfilename,experiment,year,rxNum,plothist,cutbadpix)

%mapfilename='maps/rx1_all'
%experiment='keck';
%year=2012;
%rxNum=1;
%mapfile='maps/rx1_all_comp'
%load(mapfile);
%mapfilename.map=map(:,:,:,3);
%mapfilename.x_bin=ad.t_val_deg{1};
%mapfilename.y_bin=-ad.t_val_deg{2};

[noise{1} mn{1}]=estimatenoiseinmap(mapfilename);
keyboard
if strcmp(experiment,'keck')
  if year==2012
    chflags=get_default_chflags('keck','2012');
    [p,ind]=get_array_info('20120202',[],[],[],[],chflags);
  elseif year==2013
    chflags=get_default_chflags('keck','2013');
    [p,ind]=get_array_info('20130202',[],[],[],[],chflags);
  end
  if ~strcmp(rxNum,'all')
    cutind=(p.rx==rxNum);
    p=structcut(p,cutind);
    ind=make_ind(p);
  end
elseif strcmp(experiment,'bicep2')
  chflags=get_default_chflags();
  [p,ind]=get_array_info([],[],[],[],[],chflags);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%inject noise

%ze choppa
c=0.0044;
%for uberchopper
a=8.6043e-4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noise=noise{1};
mn=mn{1};
if cutbadpix
badpix=setdiff(ind.e,ind.rgl);
noise(badpix)=NaN;
mn(badpix)=NaN;
end

if plothist
  keyboard
[h,x]=hist(noise(ind.rgl),0:5e-5:5e-3);
stairs(x,h)
hold on
plot([a a],[0 1000],'g','linewidth',2)
plot([c c],[0 1000],'r','linewidth',2)
axis([0 5e-3 0 70])
title('noise in peak normalized beam map')
legend(['rx' num2str(rxNum) 'composite bm noise'],...
    'uberchopper level','zechoppa level')
mkpng(['plots_20131217/keckcompbm_noise_hist_rx' num2str(rxNum)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [noiserms mean]=estimatenoiseinmap(mapfilename)
%find the noise rms of an input beam map
%look at ~/201300328_beamresidual/constructedmapsnoise/generateplots.m


%mapfilename='constructedmaps_v2/egaussmaps_uber.mat';

%ok, so load the beam map
if ischar(mapfilename)
load(mapfilename)

x_bin=ad.t_val{1}*180/pi;
y_bin=ad.t_val{2}*180/pi;
else
  map=mapfilename.map;
  x_bin=mapfilename.x_bin;
  y_bin=mapfilename.y_bin;  
end


%[p,ind]=get_array_info();
%p.gcp=ones(1,528);

%mask out the beam 
%mask out a 2deg square around center of map
sizeindeg=2;
maskaz= x_bin>=-sizeindeg/2 & x_bin<=sizeindeg/2;
maskel= y_bin>=-sizeindeg/2 & y_bin<=sizeindeg/2;
[mmaz mmel]=meshgrid(maskaz,maskel);
mask=mmaz&mmel;
mask=repmat(mask,[1,1,size(map,3)]);

%noisemap stack
noisemap=map;
noisemap(mask)=NaN;
%keyboard

%make noisemap into a vector
noisemap1=squeeze(reshape(noisemap,1,[],size(map,3)));

noiserms=NaN(1,size(map,3));
for ii=1:size(map,3)
  %keyboard
  %ii
  tmp=noisemap1(:,ii);
  tmp(isnan(tmp))=[];
  if ~isempty(tmp)
    [noiserms(ii) mean(ii)]=calcnoise(tmp);  
  end
end

%can i treat this as one map?? maybe??
%noisemap2=reshape(noisemap,1,[]);
%noisemap2(isnan(noisemap2))=[];
%noisrms2=calcnoise(noisemap2);
%keyboard

function [sig mu]=calcnoise(noisevector)
%fitting routine
[y,x]=hist(noisevector,100);
%stairs(x,y)
%pause
fitobject = fit(x',y','gauss1');
mu=fitobject.b1;
sig=fitobject.c1/sqrt(2);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
[y,x]=hist(noisemap1(:,1),100);  
stairs(x,y,'b')
hold on
[y,x]=hist(thing1(:,2),100);   
stairs(x,y,'r')
[y,x]=hist(thing1(:,3),100);
stairs(x,y,'g')             
[y,x]=hist(thing1(:,4),-3e-4:5e-6:3e-4);
stairs(x,y,'k')
%}

