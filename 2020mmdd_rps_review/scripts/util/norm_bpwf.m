function bpwf=norm_bpwf(bpwf)
% bpwf=norm_bpwf(bpwf)
%
% take the sum of the bpwf and return this along with the bpwf
% renormalized to unit sum

% find na, the last auto spectrum before switching
% to cross spectra:
N=length(bpwf);
na = -0.5+sqrt(2*N+0.25);

% norm to sum of 1
for j=1:length(bpwf)
  % treat TT,TE
  for i=1:2
    % take the sum over ell for each band of spec i
    norm=sum(bpwf(j).Cs_l(:,:,i));
    % store this integral
    bpwf(j).sum(:,i)=norm;
    % expand out to each ell
    norm=repmat(norm,[size(bpwf(j).Cs_l,1),1]);
    % apply to this spec
    bpwf(j).Cs_l(:,:,i)=bpwf(j).Cs_l(:,:,i)./norm;
  end
  % treat EE,BB
  for i=3:4
    % take the sum over ell for the two spec
    norm=sum(bpwf(j).Cs_l(:,:,[i,i+2]));
    % take the sum over ee->ee/bb
    norm=sum(norm,3);
    % store this integral
    bpwf(j).sum(:,i)=norm;
    % expand out to each ell
    norm=repmat(norm,[size(bpwf(j).Cs_l,1),1]);
    % apply to the two specs
    bpwf(j).Cs_l(:,:,i)=bpwf(j).Cs_l(:,:,i)./norm;
    bpwf(j).Cs_l(:,:,i+2)=bpwf(j).Cs_l(:,:,i+2)./norm;
  end
  if(j>na && size(bpwf(j).Cs_l,3)>6)
    % treat ET for cross spectra
    i=7;
    norm=sum(bpwf(j).Cs_l(:,:,i));
    bpwf(j).sum(:,i)=norm;
    norm=repmat(norm,[size(bpwf(j).Cs_l,1),1]);
    bpwf(j).Cs_l(:,:,i)=bpwf(j).Cs_l(:,:,i)./norm;
  end
    
end

return
