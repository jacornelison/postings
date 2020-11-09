function acc=coadd_ac_overrx(ac)
% function ac=coadd_ac_overrx(ac)
%
% written to coadd saved acs over receiver
% no channel/rx selection is done
% 
% example:
% load('1315/real_2012_filtp3_weight3_gs_dp1100_jack01.mat')
% ac=coadd_ac_overrx(ac)

if ~exist('ac','var')
  error('Need to input ac!');
end

%combcomap just concatenates ac, so have to loop through those
if iscell(ac)
  for kk=1:size(ac,1)
    for jj=1:size(ac,2)
      acc{kk,jj}=addrx(ac{kk,jj});
    end
  end
else
  acc=addrx(ac);
end

return

function acc=addrx(ac)

%initialize ac
acc=ac(1,:);

%add in each of the following rxs for each jackhalf
for jj=1:size(ac,2)
  for ii=2:size(ac,1)
    acc(1,jj)=addac_nan(acc(1,jj),ac(ii,jj));
  end
end

return

