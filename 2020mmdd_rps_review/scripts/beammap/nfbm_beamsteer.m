function [steerx steery]=nfbm_beamsteer(rxNum)
%function [steerx steery]=nfbm_beamsteer(rxNum)
%Function to figure out what the magnitude of the beam steer is
%in the near field.
%copied from ~clwong/nfbm/findaperture.m
%This isn't guarenteed to work. Also, it takes 
%data from Grant's fits from 2012...
%for data I've used this with look in ~clwong/nfbm
%function [sse c gof p newmean p2]=findaperture(rxNum)
%function aperturesize=findaperture(dist)
%function [steerx steery]=findaperture(rxNum)
%the beam steer is defined in
%inches or some unit of length
%and is the distance from the 
%center of the beam from the center of the aperture
%see postings
%http://bicep.caltech.edu/~spuder/keck_analysis_logbook/analysis/20120726_beamsteer/
%good luck!
%CLW 20 August 2014


%rxNum=0;
%dist= 11.1750 %gives aperture size=2.4338


%distance from aperture plane to near field mapping plane.
dist=12;  %gives aperture size=2.4411
if rxNum==0
  dist=13.7;
end

%fit for something
fitting=0;
%rxNum=0
load(['~clwong/nfbm/data/rx' int2str(rxNum) '_nfbm.mat'])
[pp ind]=get_array_info('20120130');
cutind=(pp.rx==rxNum);
pp=structcut(pp,cutind);
ind=make_ind(pp);
pp.theta=pp.theta+pp.drumangle;

load(['data/rx' int2str(rxNum) '_gauss_fits.mat']);

%this is some function grant wrote to clean the maps,
%you probably don't need it if your maps are clean
if ~exist('clean_maps','var')
  for ii=1:528
    clean_maps{ii}=beam_map_cleaner(my_maps{ii});
  end
end

if rxNum==10
  for ii=1:528
    [ro co pol ti]=gcp2det(ii-1,'K3_2');
    if ti~=0
    index=det2gcp(ro,co,pol,ti)+1;
    newmaps{index}=clean_maps{ii};
    end
  end
  clean_maps=newmaps;
  clear newmaps ro co pol ti;
end

if ~exist('x','var')
  x=0:0.25:0.25*(length(clean_maps{1})-1);
  y=0:0.25:0.25*(length(clean_maps{1}')-1);
end

xstep=x(2)-x(1);
ystep=y(2)-y(1);

Anew=NaN(528,7);
for ii=1:528
  Anew(ii,1)=A{ii}(1);
  Anew(ii,2)=(A{ii}(2)-1)*xstep+x(1);
  Anew(ii,3)=(A{ii}(3)-1)*ystep+y(1);
end

%mark for good maps
cutspixel=NaN(1,528);
my_flag=manual_beam_map_flags(['rx' int2str(rxNum)]);
for ro=1:8
  for co=1:8
    for ti=1:4
      index=det2gcp(ro,co,'A',ti)+1;
      cutspixel(index)=my_flag.a(ro,co,ti);
      index=det2gcp(ro,co,'B',ti)+1;
      cutspixel(index)=my_flag.b(ro,co,ti);
    end
  end
end
cutspixel(isnan(cutspixel))=1;
badpix=find(cutspixel);

%this figures out where the center of the aperture is
%moves the map, so the the map is centered on
%the center of the aperture
%this does it to the nearest pixel
%should be reworked to use interp2
xnew=NaN(528,length(x));
ynew=NaN(528,length(y));
for ii=1:528
  movex(ii)=dist*tand(pp.r(ii)).*sind(pp.theta(ii));
  movecol(ii)=round(movex(ii)/xstep);
  movey(ii)=dist*tand(pp.r(ii)).*cosd(pp.theta(ii));
  moverow(ii)=round(movey(ii)/ystep);
  xnew(ii,:)=x+movecol(ii)*xstep;
  ynew(ii,:)=y+moverow(ii)*ystep;
  xpos(ii)=Anew(ii,2)+movecol(ii)*xstep;
  ypos(ii)=Anew(ii,3)+moverow(ii)*ystep;
end


xall=min(min(xnew)):xstep:max(max(xnew));
yall=min(min(ynew)):ystep:max(max(ynew));

map=zeros(length(yall),length(xall),528);

for ii=1:528
  if ~isnan(xnew(ii,1)) && ~isnan(ynew(ii,1))
    xstart=int8((xnew(ii,1)-xall(1))/xstep+1);
    xend=int8(length(xall)-(xall(end)-xnew(ii,end))/xstep);
    ystart=int8((ynew(ii,1)-yall(1))/ystep+1);
    yend=int8(length(yall)-(yall(end)-ynew(ii,end))/ystep);

    map(ystart:yend,xstart:xend,ii)=clean_maps{ii}/Anew(ii,1);
  end
end

mapAll=zeros(size(map(:,:,1)));
for ii=1:528
  if cutspixel(ii)==0
    if Anew(ii,1)~=0 && ~isnan(Anew(ii,1))
    mapAll=mapAll+map(:,:,ii);
  end
  end
end

%imagesc(xall,yall,mapAll)

%finds the center of the aperture
center=normfit2d(xall,yall,mapAll)
Sigma = [center(4)^2  center(6)*center(5)*center(4); center(6)*center(5)*center(4) center(5)^2];
[V, D]=eigs(Sigma);
aperturesize=(sqrt(D(1,1))+sqrt(D(2,2)))/2

%mapnew=mapAll(moverow(10)+1:moverow(10)+56,movecol(10)+1:movecol(10)+56);
%figure(1)
%mapnew(mapnew<0)=0;
%imagesc(log10(mapnew))
%caxis([0 2.5])
%axis image
%title(['dist=' int2str(dist)])
%print('-dpng',['~/posting/dist' int2str(dist)])
%keyboard
%figure(2)
%hold on
%thing=normfit2d([1:56],[1:56],mapnew)
%ok=round(thing(2))
%slice=mapnew(:,ok);
%plot(mapnew(:,ok),'g')
%xlabel('Y bin')
%ylabel('Amplitude')
%title(['dist=' int2str(dist)])
%print('-dpng',['~/posting/dist' int2str(dist) '_slice'])
%figure(3)
%plot(log10(mapnew(:,ok)),'r')
%xlabel('Y bin')
%ylabel('Log Amplitude')
%dist= 11.1750

%center of the aperture
xcen=center(2);
ycen=center(3);

%the beam steer is defined in
%inches or some unit of length
%and is the distance from the 
%center of the beam from the center of the aperture

load(['data/rx' int2str(rxNum) '_gauss_fits.mat']) 
load(['data/rx' int2str(rxNum) '_centers_of_mass.mat']) 
An=NaN(528,7);
Rn=NaN(528,2);
for ii=1:528
  if cutspixel(ii)==0
  An(ii,:)=A{ii};
  Rn(ii,:)=R{ii};
end
end
%imagesc(clean_maps{ii})
%hold on
%plot(A{ii}(2),A{ii}(3),'x','linewidth',2)
%plot(R{ii}(2),R{ii}(1),'og','linewidth',2)
%legend('Gaussian fit','Center of mass fit')
%title('gcp240')
%print('-dpng','gaussianvscomfit_rx1')


%xpos=R{pix}(2);
%ypos=R{pix}(1);
%load(['data/rx' int2str(rxNum) '_aperture_center.mat'])

steerx=xpos-xcen;
steery=ypos-ycen;
steery=-steery;
steerx(badpix)=NaN;
steery(badpix)=NaN;
steermag=sqrt(steerx.^2+steery.^2);

%rectangle('Position',[center(2)-5 center(3)-5 10 10],'Curvature',[1,1])

xc=-pp.r.*sind(pp.theta)*pi/180;
xc=xc*150;
yc=pp.r.*cosd(pp.theta)*pi/180;
yc=yc*150;
xc=xc';
yc=yc';

%figure(3)
%clf
%colormap('jet')
%set(3,'Position',[10 675 600 600])
%scale=1;
%quiver(xc(ind.a),yc(ind.a),steerx(ind.a)*scale,steery(ind.a)*scale,0,'b')
%hold on
%quiver(xc(ind.b),yc(ind.b),steerx(ind.b)*scale,steery(ind.b)*scale,0,'r')
%legend('A pol','B pol')
%title(['Gaussian fit - Aperture center, rx' int2str(rxNum)])
%%%%%print('-dpng',['~/posting/apcengausfitdiff_rx' int2str(rxNum)])
%i have no idea what this is for
%some reorganization of data?
  ysa=NaN(8,8,4);
  ysb=NaN(8,8,4);
for ii=1:528  
  [ro co pol ti]=gcp2det(ii-1);
  if strcmp(pol,'A')
    ysa(ro,co,ti)=-steery(ii);
  elseif strcmp(pol,'B')
    ysb(ro,co,ti)=-steery(ii);
  end
  end
for ro=1:8
  for co=1:8
    for ti=1:4
      ys(ro,co,ti)=nanmean([ysa(ro,co,ti) ysb(ro,co,ti)]);
    end
  end
end

%vector=[-3.5 -2.5 -1.5 -0.5 0.5 1.5 2.5 3.5];

%if you want to fit the amount of beam
%steer in a tile to a slope, use this
%see here:
%http://bicep.caltech.edu/~spuder/keck_analysis_logbook/analysis/20120802_beamsteer/
if fitting==1
  s = fitoptions('Method','NonlinearLeastSquares')
  f = fittype('a*x+b','options',s);
  vector=[-3.5 -2.5 -1.5 -0.5 0.5 1.5 2.5 3.5];
  vector=vector*0.785;
  %vector=vector*7.85*0.1;
  %vector=12*tand(vector);
  vector2=[-2.5 -1.5 -0.5 0.5 1.5 2.5];
  vector2=vector2*0.785;
  
  for jj=1:2
    data=[];
    err=[];
    for ii=1:8
      data=[data nanmean([ysa(ii,:,jj) ysb(ii,:,jj)])];
      err=[err nanstd([ysa(ii,:,jj) ysb(ii,:,jj)])];
    end
    p(jj,:)=polyfit(vector,data,1);
    testdata(jj,:)=p(jj,1)*vector+p(jj,2);
    data2=data(2:end-1);
    p2(jj,:)=polyfit(vector2,data2,1);
    testdata2(jj,:)=p2(jj,1)*vector2+p2(jj,2);
    [c{jj} gof{jj}]=fit(vector',data',f);
    ss=data-testdata(jj,:);
    ss=ss.*ss;
    sse(jj)=sum(ss);
    datamean{jj}=data;
    errmean{jj}=err;
    newmean(jj)=sum(data)/8;
  end
  for jj=3:4
    data=[];
    err=[];
    for ii=1:8
      data=[nanmean([ysa(ii,:,jj) ysb(ii,:,jj)]) data];
      err=[nanstd([ysa(ii,:,jj) ysb(ii,:,jj)]) err];
    end
    p(jj,:)=polyfit(vector,data,1);
    testdata(jj,:)=p(jj,1)*vector+p(jj,2);
    [c{jj} gof{jj}]=fit(vector',data',f);
    data2=data(2:end-1);
    p2(jj,:)=polyfit(vector2,data2,1);
    testdata2(jj,:)=p2(jj,1)*vector2+p2(jj,2);
    ss=data-testdata(jj,:);
    ss=ss.*ss;
    sse(jj)=sum(ss);
    datamean{jj}=data;
    errmean{jj}=err;
    newmean(jj)=sum(data)/8;
  end


figure(1)
set(1,'Position',[10 675 600 600])
clf
jj=1;
subplot(2,2,jj)
plot(vector,testdata(jj,:),'r')
hold on
plot(vector2,testdata2(jj,:),'g')
errorbar(vector,datamean{jj},errmean{jj},'b')
legend(['y=' num2str(p(jj,1)) 'x+' num2str(p(jj,2))], ...
    ['y=' num2str(p2(jj,1)) 'x+' num2str(p2(jj,2))])
%legend('fit','fit 6 rows') 
%ext(-2.5,1.5,['mean=' num2str(mean(jj))])
xlabel('row distance from center of tile,mm')
ylabel('avg beam steer across a row, in')
axis([-3 3 -2 2])
grid on
title(['rx' int2str(rxNum) ',tile' int2str(jj) 'mean=' num2str(newmean(jj))])
%title(['rx' int2str(rxNum) ',tile1,y=' num2str(p(jj,1)) 'x+' ...
%      num2str(p(jj,2))])
jj=2;
subplot(2,2,jj)
plot(vector,testdata(jj,:),'r')
hold on
plot(vector2,testdata2(jj,:),'g')
errorbar(vector,datamean{jj},errmean{jj},'b')
legend(['y=' num2str(p(jj,1)) 'x+' num2str(p(jj,2))], ...
    ['y=' num2str(p2(jj,1)) 'x+' num2str(p2(jj,2))])
%legend('fit','fit 6 rows') 
%ext(-2.5,1.5,['mean=' num2str(mean(jj))])
xlabel('row distance from center of tile,mm')
ylabel('avg beam steer across a row, in')
axis([-3 3 -2 2])
grid on
title(['rx' int2str(rxNum) ',tile' int2str(jj) 'mean=' num2str(newmean(jj))])
jj=3;
subplot(2,2,4)
plot(vector,testdata(jj,:),'r')
hold on
plot(vector2,testdata2(jj,:),'g')
errorbar(vector,datamean{jj},errmean{jj},'b')
legend(['y=' num2str(p(jj,1)) 'x+' num2str(p(jj,2))], ...
    ['y=' num2str(p2(jj,1)) 'x+' num2str(p2(jj,2))])
%legend('fit','fit 6 rows') 
%ext(-2.5,1.5,['mean=' num2str(mean(jj))])
xlabel('row distance from center of tile,mm')
ylabel('avg beam steer across a row, in')
axis([-3 3 -2 2])
grid on
title(['rx' int2str(rxNum) ',tile' int2str(jj) 'mean=' num2str(newmean(jj))])
jj=4;
subplot(2,2,3)
plot(vector,testdata(jj,:),'r')
hold on
plot(vector2,testdata2(jj,:),'g')
errorbar(vector,datamean{jj},errmean{jj},'b')
legend(['y=' num2str(p(jj,1)) 'x+' num2str(p(jj,2))], ...
    ['y=' num2str(p2(jj,1)) 'x+' num2str(p2(jj,2))])
%legend('fit','fit 6 rows') 
%ext(-2.5,1.5,['mean=' num2str(mean(jj))])
xlabel('row distance from center of tile,mm')
ylabel('avg beam steer across a row, in')
axis([-3 3 -2 2])
grid on
title(['rx' int2str(rxNum) ',tile' int2str(jj) 'mean=' num2str(newmean(jj))])
%print('-dpng',['~/posting/dataerr_rx' int2str(rxNum)])
end
%{
clf
subplot(2,1,1)
errorbar(vector,datamean{1},errmean{1},'b')
hold on
errorbar(vector,datamean{2},errmean{2},'r')
errorbar(vector,datamean{3},errmean{3},'m')
errorbar(vector,datamean{4},errmean{4},'g')
%axis([
legend('tile 1','tile 2','tile 3','tile 4')
title(['data and err, rx' int2str(rxNum)])
subplot(2,1,2)
plot(vector,testdata(1,:),'b')
hold on
plot(vector,testdata(2,:),'r')
plot(vector,testdata(3,:),'m')
plot(vector,testdata(4,:),'g')
legend(['tile 1 ' num2str(p(1,1)) ' ' num2str(p(1,2))] ,...
    ['tile 2 ' num2str(p(2,1)) ' ' num2str(p(2,2))],...
    ['tile 3 ' num2str(p(3,1)) ' ' num2str(p(3,2))],...
    ['tile 4 ' num2str(p(4,1)) ' ' num2str(p(4,2))])
title(['linear fit, rx' int2str(rxNum)])
subplot(3,1,3)
xlim([-4 4])
plot(c{1},'b')
hold on
plot(c{2},'r')
plot(c{3},'m')
plot(c{4},'g')
legend(['tile 1 ' num2str(coeffvalues(c{1}))] ,...
    ['tile 2 ' num2str(coeffvalues(c{2}))],...
    ['tile 3 ' num2str(coeffvalues(c{3}))],...
    ['tile 4 ' num2str(coeffvalues(c{4}))])
title(['linear fit, rx' int2str(rxNum)])
print('-dpng',['plots/dataerr_rx' int2str(rxNum)])
%}
%ok, finding the average doesn't really work
%for ii=1:4
%avgnum(ii)=nanmean(datamean{ii})
%stdnum(ii)=nanstd(datamean{ii})
%end
%keyboard

%scatter plot y-ABoffset vs beam steer
%y A-B offset
%yAB=-An(ind.a,3)+An(ind.b,3);
%ys2=(-steery(ind.a)-steery(ind.b))/2;
%figure(2)
%clf
%colormap('jet')
%scatter(yAB,ys2)

%keyboard