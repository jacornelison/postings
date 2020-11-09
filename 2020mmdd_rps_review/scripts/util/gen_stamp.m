function s=gen_stamp()
% function s=gen_stamp()
% generate random string with, usefull for tempfile, exportet from farmit
% -->  s = '20121203_082431_tw4w'

  % avoid having the same string sequence on all nodes:
  % would be nice to use the 'shuffle' option, but that's
  % not available in all matlab versions
%    persistent crgn;
%    if (size(crgn,1)==0) 
%      ms=datevec(now);
%      crgn = RandStream('mt19937ar','Seed',ms(6)*1e7);
%    end
%    
%    stamp0=floor(36*rand(crgn,1,4));
%    stamp1(stamp0<10)='0'+stamp0(stamp0<10);
%    stamp1(stamp0>=10)='a'-10+stamp0(stamp0>=10);
%    stamp1=char(stamp1);
  
  
  % will create a tmp file in the /tmp space which is
  % /scratch on the odyssey nodes
  [dum,tmpfilename]=system_safe('mktemp');
  % since the scratch space is per node, still the same
  % tmpfilename can be created on different nodes. To
  % prevent that fetch the hostname...
  host = hostname();
  % and add it to the tmpfilename string. Adding two string
  % will sum the ascii codes of the two characters. We have
  % to bend this back into the ascii chart portion that is
  % standard letters and numbers, these code positions:
  asc=[48:57,65:90,97:122];
  % try for instance char(asc)
  % now bend it back into an allowed range using modulo
  rnstr = mod(tmpfilename(end-4:end-1) + host(end-4:end-1),length(asc))+1;
  % and map it onto characters:
  rnstr = char(asc(rnstr));
  
  s=[datestr(now,'yyyymmdd_HHMMSS') '_' rnstr];
  
  system_safe(['rm -f ',tmpfilename]);

return
