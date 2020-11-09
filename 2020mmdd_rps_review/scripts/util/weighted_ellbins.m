function r=weighted_ellbins(r,bpwf)
%
% Get the centroids of the band power window functions.
% Also output 16 and 84 percent points for making error bars.
%
% e.g.
% load final/0751/real_a_filtp3_weight3_gs_dp1102_jack0_matrix
% r=weighted_ellbins(r,bpwf)

for i=1:length(r)
  l=cvec(bpwf(i).l);
  Cs_l=bpwf(i).Cs_l;

  Cs_l(isnan(Cs_l))=0;

  r(i).lc=nan(size(Cs_l,2),size(Cs_l,3));
  r(i).ll=nan(size(Cs_l,2),size(Cs_l,3));
  r(i).lh=nan(size(Cs_l,2),size(Cs_l,3));

  for kk=1:size(Cs_l,3)
    for jj=1:size(Cs_l,2)
      % take the sum
      s=trapz(l,Cs_l(:,jj,kk));
      if(s>0) % don't barf on all zero BPWF (as EE->BB sometimes is)
        % normalize to unit sum (normally will be anyway)
        Cs_l(:,jj,kk)=Cs_l(:,jj,kk)/s;
        % weighted mean
        r(i).lc(jj,kk)=trapz(l,l.*Cs_l(:,jj,kk));
        % low side
        xl_ind=find(cumsum(Cs_l(:,jj,kk))>0.15865,1,'first');
        r(i).ll(jj,kk)=l(xl_ind);
        % high side
        xh_ind=find(cumsum(Cs_l(:,jj,kk))<0.84135,1,'last');
        r(i).lh(jj,kk)=l(xh_ind);
      end
    end
  end
end

return
