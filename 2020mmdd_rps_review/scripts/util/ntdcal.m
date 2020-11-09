% d=ntdcal(d,tempconvertfile)
%
% Apply NTD calibrations offline.  This is useful for
% time periods when online GCP calibration was not
% available.  The file TempConvert.cc from the GCP
% repository should be passed as the second argument.
function d=ntdcal(d,tcf)

cal=read_tempconvert_file(tcf);
xx=d.antenna0.hk.fast_voltage;
yy=zeros(size(xx));
for i=1:size(xx,2)
  if i>length(cal) || isempty(cal(i).val)
    continue
  end
  yy(:,i)=interp1(cal(i).val(:,1),cal(i).val(:,2),xx(:,i),'spline');
  yy(:,i)=exp(yy(:,i));
end
d.antenna0.hk.fast_temp=yy;

return

%%
function cal=read_tempconvert_file(fname)

f=fopen(fname,'rt');
ch=[];
while ~feof(f)
  ll=fgetl(f);
  if ~isempty(strfind(ll,'fast channel'))
    tok=regexp(ll,'fast channel (\d+)','tokens');
    if ~isempty(tok)
      ch=str2num(tok{end}{end});
    end
    continue
  end
  if ~isempty(ch)
    ll=strtrim(ll);
    if ll(1)=='{'
      cal(ch+1).val=[];
      ll=ll(2:end);
      while(true)
        tmp=sscanf(ll,'%f,%f,');
        if length(tmp)==2
          cal(ch+1).val=[cal(ch+1).val; reshape(tmp,1,2)];
        end
        ll=fgetl(f);
        ll=strtrim(ll);
        if feof(f) && isempty(ll)
          break;
        end
        if ll(1)=='}'
          break;
        end
      end
      ch=[];
    end
  end
end

return

