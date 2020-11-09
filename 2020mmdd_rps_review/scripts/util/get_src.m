function src=get_src(tag,mapopt)
% src=get_src(tag,mapopt)
%
% Get structure containing info on pnt src for excise from QUaD maps
%
% If mapopt arg provided look in there to see if we should double up
% the source list to accomodate field difference

if(~exist('tag','var'))
  tag=[];
end

if(isempty(tag))
  tag='cmb';
end

switch tag
  case 'cmb'

    % From Jones Thesis
    %src.name={'PKS0537-441','PKS0518-45','PKS0454-46'};
    %src.ra=[84.68,79.98,73.98];
    %src.dec=[-44.11,-45.82,-46.32];
    
    % same names into NED gives coords near identical to
    % PMNS coords except for double PMNS source where gives mid point
    
    % in approx decending order of importance as seen in QUaD maps
    
    % bright quasar - appears point
    src.name{1}='PKS0537-441';
    src.ra(1)=84.7098392;
    src.dec(1)=-44.0858150;
    % src.r is radius to mask to in poly sub - set to 0.2 for bright
    % sources which disrupt the poly sub and 0 otherwise
    src.r(1)=0.2;
    
    % resolved radio galaxy
    src.name{end+1}='PKS0518-45 (PICTOR A)';
    src.ra(end+1)=79.957045;
    src.dec(end+1)=-45.779019;
    src.r(end+1)=0.2;
    
    % This one ringed on Jones map but not listed in table
    % there is another source very nearby
    src.name{end+1}='PKS0514-459';
    src.ra(end+1)=78.9387667;
    src.dec(end+1)=-45.9453694;
    src.r(end+1)=0.2;

    % just on edge of QUaD map
    src.name{end+1}='PKS0454-46';
    src.ra(end+1)=73.9615517;
    src.dec(end+1)=-46.2663003;
    src.r(end+1)=0;
    
    src.name{end+1}='PKS0557-454';
    src.ra(end+1)=89.798068;
    src.dec(end+1)=-45.494505;
    src.r(end+1)=0.2;
    
    src.name{end+1}='PKS 0524-485';
    src.ra(end+1)=81.56570;
    src.dec(end+1)=-48.50972;
    src.r(end+1)=0;
    
    % Visible in 2006 region
    
    % Source PKS 0524-485 above is visible in overlap region
    
    % This source which was actually marginally visible in 2005 region
    % is now clearly detected
    src.name{end+1}='PMNJ0558-5029';
    src.ra(end+1)=89.5450;
    src.dec(end+1)=-50.4967;
    src.r(end+1)=0;
    
    % This one is visible in the quad map but real weak in PMN (100mJy)
    src.name{end+1}='PMNJ0549-5246';
    src.ra(end+1)=87.4250;
    src.dec(end+1)=-52.7694;
    src.r(end+1)=0;

    % if field diff map we double up the point sources by adding a second
    % copy shifted to smaller RA by 7.5 deg - this accounts for the fact
    % that sources in the overlap region strip down the middle can appear
    % twice - once as a positive source at the +ve RA end of the diff
    % field and once as a negative source at the negative end.
    
    if(exist('mapopt','var'))
      
      % this works for real maps
      if(mapopt.fd==1)
	src.ra=[src.ra,src.ra-7.5];
	src.dec=[src.dec,src.dec];
	src.name=[src.name,src.name];
	src.r=[src.r,src.r];
      end

      % for sim field diff happens in makesim so we need to look there
      if(isfield(mapopt,'simopt'))
	if(mapopt.simopt{1}.fd==1)
	  src.ra=[src.ra,src.ra-7.5];
	  src.dec=[src.dec,src.dec];
	  src.name=[src.name,src.name];
	  src.r=[src.r,src.r];
	end
      end
      
    end
    
  case 'rcw38'
    src.name{1}='RCW38';
    src.ra(1)=134.77;
    src.dec(1)=-47.51;
    src.r=0.5;

    % extra sources in this field
    src.name{end+1}='RCW38a';
    src.ra(end+1)=135.4721;
    src.dec(end+1)=-47.7292;
    src.r(end+1)=0.2;
    
    src.name{end+1}='RCW38b';
    src.ra(end+1)=135.5609;
    src.dec(end+1)=-48.6892;
    src.r(end+1)=0.2;
    
    src.name{end+1}='RCW38c';
    src.ra(end+1)=135.8570;
    src.dec(end+1)=-48.4492;
    src.r(end+1)=0.2;
  case 'mat6a'
    src.name{1}='mat6a';
    src.ra(1)=12.16767*15;
    src.dec(1)=-62.846;
    src.r(1)=0.2;

    % extra source in this field
    %src.name{end+1}='mat6a';
    %src.ra(end+1)=12.16767*15;
    %src.dec(end+1)=-62.846;
    %src.r(end+1)=0.2;    
  case 'g326.3'
    src.name{1}='g326.3';
    src.ra(1)=15.8833*15;
    src.dec(1)=-56.1667;
    src.r(1)=0.5;
  case 'carneb'
    src.name{1}='carneb';
    src.ra(1)=10.73*15;
    src.dec(1)=-59.65;
    src.r(1)=1.0;
  case 'galsurv'
    % use src_coords to get field ra
    % scans sets are 0.5deg in dec at constant az range, so we don't
    % want a continuous line in RA and dec
    % things get annoying because we work 'up' the field in dec
    dec=-62:0.02:-20;
    dec=fliplr(dec);
    for i=1:numel(dec)
      % limit the dec range of each mask to 0.5 degrees
      if(dec(i)-floor(dec(i))==0|dec(i)-floor(dec(i))==0.5)
	% we have bright and faint arms...
	[racen(i),deccen(i)]=src_coords(sprintf('galp%2.2f',dec(i)));
	[racenf(i),deccenf(i)]=src_coords(sprintf('galpf%2.2f',dec(i)));
	rafix=racen(i);
	rafixf=racenf(i);
      else
	racen(i)=rafix;
	deccen(i)=dec(i);
	racenf(i)=rafixf;
	deccenf(i)=dec(i);
      end
      
      src.name{i}=sprintf('galp%2.2f',deccen(i));
    end

    src.ra=[racen racenf];
    src.dec=[deccen deccenf];
    % hard-wire in 1.5 degrees on sky for mask range
    src.r=1.5*ones(size(src.ra)).*cosd(src.dec);

    % rcw38 is in galsurv field
    src.name{end+1}='RCW38';
    src.ra(end+1)=134.77;
    src.dec(end+1)=-47.51;
    src.r(end+1)=0.5;

    % extra sources in this field
    src.name{end+1}='RCW38a';
    src.ra(end+1)=135.4721;
    src.dec(end+1)=-47.7292;
    src.r(end+1)=0.2;
    
    src.name{end+1}='RCW38b';
    src.ra(end+1)=135.5609;
    src.dec(end+1)=-48.6892;
    src.r(end+1)=0.2;
    
    src.name{end+1}='RCW38c';
    src.ra(end+1)=135.8570;
    src.dec(end+1)=-48.4492;
    src.r(end+1)=0.2;

    % sources derived from basic watershed extraction
    % no sources found outside 1.5 degrees from galactic plane in low
    % el part of obs
    
    % first eleven from 070319_filtp1_weight0_jack0.mat
    src.name{end+1}='galsrc1';
    src.ra(end+1)=259.9839;
    src.dec(end+1)=-38.9043;
    src.r(end+1)=0.8903;
    
    src.name{end+1}='galsrc2';
    src.ra(end+1)=253.5895;
    src.dec(end+1)=-45.2476;
    src.r(end+1)=0.5015;
    
    src.name{end+1}='galsrc3';
    src.ra(end+1)=254.8679;
    src.dec(end+1)=-40.1481;
    src.r(end+1)=0.4918;
    
    src.name{end+1}='galsrc4';
    src.ra(end+1)=257.3767;
    src.dec(end+1)=-41.5785;
    src.r(end+1)=0.3868;
    
    src.name{end+1}='galsrc5';
    src.ra(end+1)=262.2926;
    src.dec(end+1)=-36.6085;
    src.r(end+1)=0.5140;
    
    src.name{end+1}='galsrc6';
    src.ra(end+1)=255.1454;
    src.dec(end+1)=-40.5186;
    src.r(end+1)=0.3192;
    
    src.name{end+1}='galsrc7';
    src.ra(end+1)=243.7678;
    src.dec(end+1)=-49.8168;
    src.r(end+1)=0.2523;
    
    src.name{end+1}='galsrc8';
    src.ra(end+1)=250.0125;
    src.dec(end+1)=-48.8511;
    src.r(end+1)=0.2646;
    
    src.name{end+1}='galsrc9';
    src.ra(end+1)=259.5982;
    src.dec(end+1)=-39.2811;
    src.r(end+1)=0.2459;
    
    src.name{end+1}='galsrc10';
    src.ra(end+1)=245.6583;
    src.dec(end+1)=-50.7190;
    src.r(end+1)=0.2257;
    
    src.name{end+1}='galsrc11';
    src.ra(end+1)=254.4544;
    src.dec(end+1)=-40.3530;
    src.r(end+1)=0.2646;

    % next seven from 070322_filtp1_weight0_jack0.mat
    src.name{end+1}='galsrc12';
    src.ra(end+1)=134.7850;
    src.dec(end+1)=-47.4192;
    src.r(end+1)=0.8649;
    
    src.name{end+1}='galsrc13';
    src.ra(end+1)=134.7836;
    src.dec(end+1)=-47.5550;
    src.r(end+1)=0.4787;
    
    src.name{end+1}='galsrc14';
    src.ra(end+1)=135.5783;
    src.dec(end+1)=-48.6643;
    src.r(end+1)=0.2985;
    
    src.name{end+1}='galsrc15';
    src.ra(end+1)=134.8655;
    src.dec(end+1)=-43.7306;
    src.r(end+1)=0.3141;
    
    src.name{end+1}='galsrc16';
    src.ra(end+1)=135.8717;
    src.dec(end+1)=-48.4236;
    src.r(end+1)=0.2394;

    src.name{end+1}='galsrc17';
    src.ra(end+1)=135.4806;
    src.dec(end+1)=-47.7131;
    src.r(end+1)=0.2185;
    
    src.name{end+1}='galsrc18';
    src.ra(end+1)=141.3923;
    src.dec(end+1)=-47.4963;
    src.r(end+1)=0.1784;
    
end

return
