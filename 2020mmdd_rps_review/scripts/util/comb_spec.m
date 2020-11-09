function ro=comb_spec(r1,r2)
% simple combine of two sets of spectra to single set
% 
% e.g. to continue to allow specjack to work
% load final/0751x1651/real_a_filtp3_weight3_gs_dp1102_jack0_real_a_filtp3_weight3_gs_dp1000_jack0_pureB_matrix
% ro(1)=r(1); ro(1).w1=[]; ro(1).w2=[];
% ro(2)=comb_spec(r(2),r(3));
% ro(3)=comb_spec(r(4),r(5));

for i=1:size(r1.real,2) % for each spectrum
  for j=1:length(r1.l) % for each bandpower
    
    % get the sim bandpowers
    s(:,1)=r1.sim(j,i,:);
    s(:,2)=r2.sim(j,i,:);
    
    % take their covariance
    c=cov(s);
    
    % take the weight as column sum of inv cov mat
    w=sum(inv(c));
 
    % normalize to unit sum
    w=w./sum(w);
    
    % craziness supression recipe
    % taken from reduc_final_comb on B1 branch
    if(1)
      while(any(w<0))
        c=c+diag(diag(c)*1e-2);    
        w=sum(inv(c));
        w=w./sum(w);
      end
    end
    
    w1(j,i)=w(1);
    w2(j,i)=w(2);
  end
end

% apply the weights
ro.l=r1.l;
ro.expv=r1.expv.*w1+r2.expv.*w2;
ro.sim=r1.sim.*repmat(w1,[1,1,size(r1.sim,3)])+r2.sim.*repmat(w2,[1,1,size(r2.sim,3)]);
if(isfield(r1,'simr'))
  ro.simr=r1.simr.*repmat(w1,[1,1,size(r1.simr,3)])+r2.simr.*repmat(w2,[1,1,size(r2.simr,3)]);
end
ro.sig=r1.sig.*repmat(w1,[1,1,size(r1.sig,3)])+r2.sig.*repmat(w2,[1,1,size(r2.sig,3)]);
ro.noi=r1.noi.*repmat(w1,[1,1,size(r1.noi,3)])+r2.noi.*repmat(w2,[1,1,size(r2.noi,3)]);
ro.db=r1.db.*w1+r2.db.*w2;
ro.real=r1.real.*w1+r2.real.*w2;
ro.rwf=r1.rwf.*w1+r2.rwf.*w2;

% store the weights for reference
ro.w1=w1;
ro.w2=w2;

return
