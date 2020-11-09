function create_scancheck(tags,location,update)
% create_scancheck(tag)
%
% Creates scancheck_v1-v3 in the file location specified
%
% Example:
% create_scancheck('B20070706_cmb16-2C_B_135')
%
% Options:
% tags: tags you want to analyze.  format the same as
%   from get_tags.
% location: for any other location than '/data/tods/'.  
% update: this program will check to see if scancheck_v*s 
%   already exist in the directory. if update==0 (default) 
%   it wont write over the existing files.  if update==1 then
%   it will write over them. 
%
% SAS 12/8/09

if(~exist('tags','var'))
  tags={};
end

if(isempty(tags))
  tags={'B20060515_cmb10_1A_B_315'};
end

if(~isa(tags,'cell'))
  tags={tags};
end

if(~exist('location','var'))
  location=[];
end

if(isempty(location))
  location='/data/tods';
end

if(~exist('update','var'))
  update=0;
end

n_tags=size(tags,2);

for k=1:n_tags
   tag=tags{k};

   [p,ind]=get_array_info(tag);
   date=tag(2:9);
   phase=tag(end-4:end);
   cmb=tag(11:15);
   type=tag(17:18);

   if ~strcmp(location(end:end),'/') 
    location=strcat(location,'/');
   end
   dir_tod=sprintf('%s%s_%s-%s/%s/',location,date,cmb,type,phase);

   sc_existance=zeros(3);
   %check to see if scancheck_v*s exist
   if(exist(strcat(dir_tod,'scancheck_v1')))
      sc_existance(1)=1;
   end
   if(exist(strcat(dir_tod,'scancheck_v2')))
      sc_existance(2)=1;
   end
   if(exist(strcat(dir_tod,'scancheck_v3')))
      sc_existance(3)=1;
   end
 
   %read in scancheck and make scancheck_v*s
   scancheck=readdir('scancheck',dir_tod);
   scancheck=uint32(scancheck);
   scancheck_glitch=scancheck*0;
   scancheck_weather=scancheck*0;
   scancheck_dead=scancheck*0;
   n_hscans=size(scancheck,2);

   % read in scanlist.txt 
   fid=fopen(strcat(dir_tod,'scanlist.txt'));
   result=textscan(fid,'%n %n %n %n %n','commentStyle','#');
   fclose(fid);
   fields={'set', 'num', 'start', 'stop', 'inc'};
   fs=cell2struct(result, fields, 2);

   scind=ind.gl150(find(ind.gl150 ~= 82)); %removes a bad channel

   % will base the making on nestled ifs so that we have consistant
   % scanchecks (i.e. we dont make _v1 without updating or making _v2
   % and _v3.)
  if sc_existance(3)==0 | update==1
    if sc_existance(2)==0 | update==1
      if sc_existance(1)==0 | update==1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     disp('Calculating Scancheck_v1');
      % glitch removal
      %read in utcfast
      utcfast=readdir('antenna0.time.utcfast',dir_tod,1);
      n_slow=size(utcfast,1);
      start=mjd2int(utcfast(1));
      stop=mjd2int(utcfast(n_slow));
 
      out=zeros(144);  %number of bad half scans
      hit_count=zeros(144); %number of cosmic hits
      
      for i=1:size(ind.l,2)/2;   %over all light channels
        boloname1=sprintf('antenna0.bolo.mag%d',ind.l(2*i-1));
        boloname2=sprintf('antenna0.bolo.mag%d',ind.l(2*i));
        bolo1=readInterval(boloname1,start,stop,0,0,0,dir_tod);
        bolo2=readInterval(boloname2,start,stop,0,0,0,dir_tod);
        n=size(bolo1,1);
        t=(0.0:n)/50.0;
        q=find(bitget(scancheck(:,ind.l(2*i-1)),2)); %this is the skew/kurtosis
        nq=sum(bitget(scancheck(:,ind.l(2*i-1)),2));
        out(ind.l(2*i-1):ind.l(2*i))=100.*nq/n_hscans;

        for j=1:n_hscans   %over all halfscans
	  startq=fs.start(j)*5.0;  %*5 to reference it to fast tod
	  stopq=fs.stop(j)*5.0;
	  bolo1q=bolo1(startq:stopq);
	  bolo2q=bolo2(startq:stopq);
	  tq=t(startq:stopq);
	  
	  thresh1=-7*std(bolo1q-smooth(bolo1q,7));
          tmp1=find((bolo1q-smooth(bolo1q,7))<thresh1);
	  if size(tmp1,2)>0
	    hit_count(ind.l(2*i-1))=hit_count(ind.l(2*i-1))+1;
            scancheck_glitch(ind.l(2*i-1),j)=1;
	  end

	  thresh2=-7*std(bolo2q-smooth(bolo2q,7));
          tmp2=find((bolo2q-smooth(bolo2q,7))<thresh2);
	  if size(tmp2,2)>0
	    hit_count(ind.l(2*i))=hit_count(ind.l(2*i))+1;
            scancheck_glitch(ind.l(2*i),j)=1;
	  end
        end
      end 
      
      disp('Writting Scancheck_v1 to file');
      scancheck_v1=scancheck+4*scancheck_glitch;
      write_dir(scancheck_v1,'scancheck_v1',dir_tod);

      else
      disp('Skipping Scancheck_v1');
      scancheck=readdir('scancheck_v1',dir_tod);
	
    end  % for existance of _v1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Calculating Scancheck_v2');
      %weather cut
    for k=1:max(fs.set)
      start_ss=min(find(fs.set==k)); %the scanset start index
      stop_ss=max(find(fs.set==k));  %the scanset stop index
      weathercut=sum(bitget(scancheck(scind,start_ss),1));
      if weathercut>=4 
        scancheck_weather(:,start_ss:stop_ss)=1;
      end
    end

    disp('Writting Scancheck_v2 to file');
    scancheck_v2=scancheck_v1+8*scancheck_weather;
    write_dir(scancheck_v1,'scancheck_v2',dir_tod);

   else

    disp('Skipping Scancheck_v2 and v1');
    scancheck=readdir('scancheck_v2',dir_tod);
   end  %for existance of _v2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp('Calculating Scancheck_v3');
%dead channel cut definitions
   deadchannels=[111,104,88,128];
   n_deadchannels=size(deadchannels,2);

   beginning_date=['20080810','20071202','20080808','20070918'];
   beginning_phase=['B','B','B','E'];
   ending_date=['20081109','20080115','20080808','20080115'];
   ending_phase=['C','F','F','F'];

   for i=1:n_deadchannels
    if date > beginning_date(i) and date < ending_date(i) 
      scancheck_dead(dead_channels(i),:)=1;
    end
    if date==beginning_date(i) and phase>=beginning_phase(i)
      scancheck_dead(dead_channels(i),:)=1;
    end 
    if date==ending_date(i) and phase<=ending_phase(i)
       scancheck_dead(dead_channels(i),:)=1;
    end
   end
    
   disp('Writting Scancheck_v3 to file');
   scancheck_v3=scancheck_v2+16*scancheck_dead;
   write_dir(scancheck_v3,'scancheck_v3',dir_tod);

   else
     disp('Skipping Scancheck_v3, v2 and v1. Nothing to do.')

 end %for exitance of _v3



end
return
