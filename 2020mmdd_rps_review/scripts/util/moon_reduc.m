function moon_reduc(t1,t2,dks)
% Moon reduc
% t1 = start time
% t2 = stop time
% dks = nominal dk angle, to nearest degree
% full example:
% moon_reduc('2010-jan-01:02:03:04','2010-jan-01:05:06:07',90)

% *** There is still some user specific shit buried in here - you'll have
% to edit the code to choose where you want the data products saved.  Data
% products include:
% choose_a_directory/moon_*.mat, named according to the time stamp
% choose_a_directory/beams_bicep2_obs_moon_*.csv, named according to the
% time stamp

% There is also a hard-coded boresight offset coded in here.  In the future
% that will be accounted for in the full pointing model.

% Ver 1.0 RWA 20101216

files=list_arc_files('arc/',t1,t2);

% load the data
% apply pointing model
% Extract the structure elements of interest
% Cut the data within the structure

parfor fl=1:length(files)
    fl
    load_file=char(strcat('arc/', files(fl)));
    % Load the data
    b=load_arc(load_file);
    b=apply_pointing(b);
    b=of_interest(b,load_file);
    d(fl)=cut_struc(b,1:2:length(b.fb(:,528)));
end 
  
% Concatenate data from all of the different files:
fracUTC=vertcat(d(:).fracUTC);
ts=vertcat(d(:).fb);
source_az=vertcat(d(:).source_az);
source_el=vertcat(d(:).source_el);

% Grab nominal pointing centers
[p,k]=ParameterRead('aux_data/beams/beams_bicep2_obs_20100709.csv');

% A boresight offset may be necessary
% bor_off=[-0.7 1.38];
bor_off=[-0.58 1.38];  % Second attempt after looking at data 12/13/2010

% Concatenate position information (d.pointing.hor.az, ...) over
% different files:
az=vertcat(d(:).az);
el=vertcat(d(:).el);
dk=vertcat(d(:).dk);

% cut on dk:
cc=round(dk)==dks;
ts=ts(cc,:);
el=el(cc);
az=az(cc);
source_el=source_el(cc);
source_az=source_az(cc);

% At every moment of time, find the range and azimuth from the
% boresight center to the source:
[bor_r bor_theta]=distance(el,az,source_el,source_az);

% Construct latitude and longitudinal coordinates of the telescope boresight on the sphere
% centered on the source for every moment in time. 
[lat lon]=reckon(0,0,bor_r,bor_theta);

% Construct pixel coordinates on that same sphere, with a possible
% offset 
[pxel pxaz]=reckon(bor_off(2),bor_off(1),p.r,p.theta-90+dks);
pxaz(pxaz>90)=pxaz(pxaz>90)-360;
pxaz=circshift(pxaz,-1);
pxel=circshift(pxel,-1);

% Specify resolution of the map:
el_res=.0625;
elsamp=floor(.9*(max(lat)-min(lat))/el_res); 
azsamp=floor(.5*(max(lon)-min(lon))/el_res);

% Construct map bins:
[n x_bin]=hist(lon,azsamp);
[n y_bin]=hist(lat,elsamp);

% Make the map
for ii=1:528
    map(:,:,ii)=grid_map(lon, lat, ts(:,ii),x_bin,y_bin);
end

% specify window half-width:
win_hw=10;

% Construct windowed map bins:
x_win=-abs(x_bin(2)-x_bin(1))*win_hw:abs(x_bin(2)-x_bin(1)):abs(x_bin(2)-x_bin(1))*win_hw;
y_win=-abs(y_bin(2)-y_bin(1))*win_hw:abs(y_bin(2)-y_bin(1)):abs(y_bin(2)-y_bin(1))*win_hw;

% Find indices cooresponding to the closest approach to pixel centroid:
parfor ii=1:528
    [mm xc_ind(ii)]=min(abs(x_bin-pxaz(ii)));
    [mm yc_ind(ii)]=min(abs(y_bin-pxel(ii)));
end

winm=NaN(length(x_win),length(y_win),528);
A=zeros(7,528);

% Construct and fit windowed maps:
n=length(winm(:,1,1));
for ii=1:528
    if (yc_ind(ii) > win_hw) && (xc_ind(ii) > win_hw) && (yc_ind(ii)<(elsamp-win_hw)) && (xc_ind(ii)<(azsamp-win_hw))
        winm(:,:,ii)=map(yc_ind(ii)-win_hw:yc_ind(ii)+win_hw,xc_ind(ii)-win_hw:xc_ind(ii)+win_hw,ii);
        winm(:,:,ii)=nan_interp(winm(:,:,ii));
        sb=mean(winm(1:3,:,ii),1);
        sb2=mean(winm(:,1:3,ii),2);
        pf=polyfit(1:n,sb,2);
        pf2=polyfit((1:n)',sb2,2);
        pf(3)=0;
        subp=polyval(pf,(1:n));
        subp2=polyval(pf2,(1:n)');
        pln=zeros(n,n);
        pln = pln + repmat(subp,n,1) + repmat(subp2,1,n);
        winm(:,:,ii)=winm(:,:,ii)-pln;
        [A(:,ii) Z]=normfit2d(x_win, y_win, winm(:,:,ii));
        pk(ii)=max(max(Z));
    else
        disp('crapout')
        A(:,ii)=NaN;
    end 
end

% Construct pixel r's and theta's for saving to csv file:
pf=p;
for ii=1:528
    gcp=p.gcp(ii);
    if gcp==0
        gcp=528;
    end
    [pf.r(ii),pf.theta(ii)]=distance(pxel(gcp)+A(3,gcp)-bor_off(2),pxaz(gcp)+A(2,gcp)-bor_off(1),0,0);
    pf.theta(ii)=pf.theta(ii)-90-dks;
end

% make sure they're within range:
pf.theta(pf.theta>360)=pf.theta(pf.theta>360)-360;
pf.theta(pf.theta<0)=pf.theta(pf.theta<0)+360;

% write to disk:
ParameterWrite(strcat('pxfits/beams_bicep2_obs_moon_',files{1}(1,1:15),'.csv'),pf,k);

save(['beammaps/moon/201012/moon' files{1}(1,1:15)],'map','x_bin','y_bin','winm','A','x_win','y_win')
return

function b=apply_pointing(b)
lat=-89.99106667;
lon=-44.65;

samprate=length(b.antenna0.pmac.fast_az_pos)/size(b.antenna0.tracker.horiz_mount,1);
% Generate pointing model
b.t=b.antenna0.time.utcfast(:,1);
b.datenum=utc2datenum(b.antenna0.time.utcfast);
b=arcvar_to_azel_ffflat(b);
pm=get_pointing_model(mean(b.t),b);
b=invpointing_model(b,pm,1);
[b.pointing.cel.ra,b.pointing.cel.dec]=azel2radec(...
    b.pointing.hor.az,b.pointing.hor.el,b.t,lat,lon);
b.pointing.cel.dk=b.pointing.hor.dk-parallactic(...
    b.pointing.hor.az,b.pointing.hor.el,b.pointing.cel.dec,lat);
b.az=b.pointing.hor.az;
b.el=b.pointing.hor.el;
b.dk=b.pointing.hor.dk;

% Get source ra and dec.  Equatorial geocentric apparent position:
source_ra=repmat(double(b.antenna0.tracker.equat_geoc(:,1))'/3.6E6,samprate,1);
b.source_ra=cvec(source_ra);
source_dec=repmat(double(b.antenna0.tracker.equat_geoc(:,2))'/3.6E6,samprate,1);
b.source_dec=cvec(source_dec);

% Get source az and el.  Horizontal geocentric apparent position:
source_az=repmat(double(b.antenna0.tracker.horiz_geoc(:,1))'/3.6E6,samprate,1);
b.source_az=cvec(source_az);
source_el=repmat(double(b.antenna0.tracker.horiz_geoc(:,2))'/3.6E6,samprate,1);
b.source_el=cvec(source_el);
source_pa=repmat(double(b.antenna0.tracker.horiz_geoc(:,3))'/3.6E6,samprate,1);
b.source_pa=cvec(source_pa);

%[b.source_ra,b.source_dec]=azel2radec(b.source_az,b.source_el,b.t,lat,lon);

function d=of_interest(b,load_file)

% Assign var names
b.fb=double(b.mce0.data.fb);

% Usual trick: Shift gcp indeces from 1:528 to 0:527, then make 0=528
b.fb=circshift(b.fb,[0,-1]);

% Generate cut:
inc=interp1(1:20:length(b.fb(:,1)),double(b.array.frame.features),1:length(b.fb(:,1)))';
cc=inc==3 & b.mce0.header.clock_counter~=0;

b.fracUTC=(double(b.antenna0.time.utcfast(:,1))+double(b.antenna0.time.utcfast(:,2))/86400);
b.datenum=utc2datenum(b.antenna0.time.utcfast);
% Cut used variables:
[d isem]=cut_struc(b,cc);

if isem==1
    disp([load_file ' is empty'])
    return
end


function [d isem]=cut_struc(b,cc)

d.fb=b.fb(cc,:);
d.eloff=b.eloff(cc);
d.azoff=b.azoff(cc);
d.pointing.cel.ra=b.pointing.cel.ra(cc);
d.pointing.cel.dec=b.pointing.cel.dec(cc);
d.az=b.az(cc);
d.el=b.el(cc);
d.dk=b.dk(cc);
d.fracUTC=b.fracUTC(cc);
d.source_ra=b.source_ra(cc);
d.source_dec=b.source_dec(cc);
d.source_az=b.source_az(cc);
d.source_el=b.source_el(cc);
d.source_pa=b.source_pa(cc);

if isempty(d.fb(:,528))
    disp('empty file')
    isem=1; 
    return
else
    isem=0;
end

function [px win]=sortstruc(map,winm,A,az_bin,el_bin,pk,xm,ym)
r1=.5; 
[AZ EL]=meshgrid(az_bin,el_bin);
% Input data into old data structure format for plotting later:
for ii=1:528
    if ii==528, ind=0; else ind=ii; end
    if ~isnan(A(:,ii))
        [ro co pol ti]=gcp2det(ind);
        WI=(AZ-A(2,ii)).^2 + (EL-A(3,ii)).^2;
        WIN=WI<r1^2;
        pow=sum(sum(map(WIN)));
        varxy=A(6,ii)*A(5,ii)*A(4,ii);
        Sigma= [A(4,ii).^2 varxy ; varxy A(5,ii).^2];
        [V, D]=eigs(Sigma);
        sig=(sqrt(D(1,1))+sqrt(D(2,2)))/2; 
        gint=D(1,1)*D(2,2)*pk(ii);
        if strcmp(pol,'A')
            px(ro,co,ti).a_SIG=[sqrt(D(1,1)),sqrt(D(2,2))];
            px(ro,co,ti).a_sig=sig;
            px(ro,co,ti).a_e=(sqrt(D(1,1))-sqrt(D(2,2)))/sig;
            px(ro,co,ti).a_cent=[A(2,ii)+xc(ii),A(3,ii)+yc(ii)];
            px(ro,co,ti).a_vec=V;
            px(ro,co,ti).a_pk=A(7,ii);
            px(ro,co,ti).a_pk2=A(1,ii);
            px(ro,co,ti).a_pow=pow;
            px(ro,co,ti).a_map=map(:,:,ii);
            px(ro,co,ti).a_gint=gint;
            win(ro,co,ti).a_map=winm(:,:,ii);
            
        elseif strcmp(pol,'B')
            px(ro,co,ti).b_SIG=[sqrt(D(1,1)),sqrt(D(2,2))];
            px(ro,co,ti).b_sig=sig;
            px(ro,co,ti).b_e=(sqrt(D(1,1))-sqrt(D(2,2)))/sig;
            px(ro,co,ti).b_cent=[A(2,ii)+xc(ii),A(3,ii)+yc(ii)];
            px(ro,co,ti).b_vec=V;
            px(ro,co,ti).b_pk=A(7,ii);
            px(ro,co,ti).b_pk2=A(1,ii);
            px(ro,co,ti).b_pow=pow;
            px(ro,co,ti).b_map=map(:,:,ii);
            px(ro,co,ti).b_gint=gint;
            win(ro,co,ti).b_map=winm(:,:,ii);
        end 
        
    else
        disp('Bad Fit')
    end
end

for ro=1:8
    for co=1:8
        for ti=1:4
            if ~isempty(px(ro,co,ti).a_cent) && ~isempty(px(ro,co,ti).b_cent)
                sig=mean(px(ro,co,ti).a_sig(:));
                px(ro,co,ti).dif_e=abs((px(ro,co,ti).b_e-px(ro,co,ti).a_e)/2);
                px(ro,co,ti).dif_p=(((px(ro,co,ti).a_cent(1)-px(ro,co,ti).b_cent(1))^2+(px(ro,co,ti).a_cent(2)-px(ro,co,ti).b_cent(2))^2)^(1/2))/sig;
                px(ro,co,ti).dif_fwhm=(mean(px(ro,co,ti).a_sig(:))-mean(px(ro,co,ti).b_sig(:)))/(mean(px(ro,co,ti).a_sig(:))+mean(px(ro,co,ti).b_sig(:)))/2;
                px(ro,co,ti).dif_pd=(((px(ro,co,ti).a_cent(1)-px(ro,co,ti).b_cent(1))^2+(px(ro,co,ti).a_cent(2)-px(ro,co,ti).b_cent(2))^2)^(1/2));
            end
        end
    end
end


% Filtering shit:
% fts=ts;
% parfor ii=1:528
%     ii
%     pf(:,ii)=polyfit((1:length(ts))',ts(:,ii),3);
%     subp=polyval(pf(:,ii),(1:length(ts))');
%     fts(:,ii)=ts(:,ii)-subp;
% end

% [bf af]=butter(1, 1/4000, 'high');
% fts=filtfilt(bf,af,ts(:,343));
