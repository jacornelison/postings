function xind = aps_getxind(aps,ind_a1, ind_a2)
% function xind = aps_getxind(aps,ind_a1, ind_a2)
% get the index of the cross spectrum in aps structure
% where ind_a1 and ind_a2 mark the postion of the two 
% corresponding auto spectra. N is the total number of

if(isstruct(aps))
  N = size(aps,1);
else
  N=aps;
end

na = -0.5+sqrt(2*N+0.25);
nc_c=na+1;

if(ind_a1>na | ind_a2>na)
  error('indecies out of range');
end

xind=[];
if na>1
  for j=1:na-1
    for c=j+1:na
      if j==ind_a1 && c==ind_a2
        xind=nc_c;
        break;
      end
      nc_c=nc_c+1;
    end
  end
end

if(isempty(xind))
  xind=ind_a1;
end

return