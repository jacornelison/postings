function [supercon, unlock, xtalk, nodata, dropped, antenna] = latching(en,endat,lcdat,p,ind)

badgain=find(abs(en.g(1,ind.gl,2))<5e2 | en.g(1,ind.gl,3)>1e4);

%Let's figure out what all of the bad gains are coming from!

%Classifications:
% 1. Superconducting... huge variance, stable line for the fist little while
% 2. Unlocked.... also huge std
% 3. Just crosstalk.  very low ampitude.
% 4. No data.  In one of the data masks, but not in fp_data
% 5. Superconducting because PLC goes too low.  good at first.  bad at end.
% 6. Antenna problems.  good in plc, no amplitude in elnod.

supercon=[];
unlock=[];
xtalk=[];
nodata=[];
dropped=[];
antenna=[];

for i=1:length(badgain)
  pix=ind.gl(badgain(i));
  enstd=std(endat(1,pix).y);
  if length(lcdat(1,pix).idet)==0
    if enstd == 0
       nodata = [nodata,pix];  % no data
    elseif enstd<10
       xtalk = [xtalk, pix];   %very little response.  usually means xtalk
    end
    continue
  end

  len=min(100,length(lcdat(1,pix).idet));
  plcstdb=std(lcdat(1,pix).idet(1:len));   %behaves nicely when normal ?
  plcstde=std(lcdat(1,pix).idet(end-len+1:end));  %behaves when supercon?
  p=polyfit(1:len,lcdat(1,pix).idet(1:len)',1);  %linear when normal
  %does it fit?
  fit=std(lcdat(1,pix).idet(1:len)'-polyval(p,1:len));
  enstd=std(endat(1,pix).y);

  if plcstdb < 1e-8     %little response
    xtalk=[xtalk pix];
  elseif isnan(plcstdb)
    if enstd == 0
      nodata=[nodata pix];
    elseif enstd < 10
      xtalk = [xtalk pix];
    end
  else

  if fit < 1e-7 %linear in first 100 (when "normal")
    if p(2) > 1e-4
      supercon = [supercon pix];    % has the slope of supercon in the first bit
    else
      if plcstde > 1e-6
         dropped = [dropped pix];   % normal, then crazy
      else
         antenna = [antenna pix];   % looks good in plc
      end
    end
   else
     
   %taken out everything else...
   unlock = [unlock, pix];    % not linear in the beginning of plc.
   end
   end

end

%Make plots of all of the bad elnod response channels
if(0)
for i=1:length(badgain)
  
  indx=ind.gl(badgain(i));
  [r c]=ind2rc(indx);
  
  clf;
  plot(lcdat(1,indx).idet);
  title(['Partial Load Curve - gcp' num2str(indx) ' r' num2str(r) 'c' num2str(c)]);
  xlabel('Sample')
  ylabel('idet')
  saveas(gcf,[folder 'lcdat_' num2str(indx) '.png'])

  clf;
  plot(endat(1,indx).y);
  title(['Elnod - gcp' num2str(indx) ' r' num2str(r) 'c' num2str(c)]);
  xlabel('Sample')
  ylabel('FB')
  saveas(gcf,[folder 'endat_' num2str(indx) '.png'])

end
end

end
