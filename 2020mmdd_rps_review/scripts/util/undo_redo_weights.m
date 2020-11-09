function ac=undo_redo_weights(subac,subhs,realhs,ind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ac=undo_redo_weights(subac,subhs,realhs,ind,acsparse,acpack) takes the flavor
% substitution ac structure, subac, divides  out the mean weight and variance across
% halfscans from the flavor substitution's mapopt.hs, subhs, and multiplies by the mean
% weight and variance across halfscans from the real pairmap mapopt.hs.  ind is the
% index structure returned from get_array_info from the flavor substitution tag.  

% deprojection templates should be stored in a cell array in ac structure
% as opposed to individual field names for each template.  As of 07/02/13, this is
% automatically taken care of in reduc_coaddpairmas before undo_red_weights is called.

% calculate the mean weights and variances from the substituted map and real map
% for pair sum and pair diff. Pair sum channel pairs are stored in ind.a, pair diff
% channels in ind.b
sumsubw=mean(subhs.w(:,ind.a),1); sumsubv=mean(subhs.s(:,ind.a).^2,1);
sumrealw=mean(realhs.w(:,ind.a),1); sumrealv=mean(realhs.s(:,ind.a).^2,1);
difsubw=mean(subhs.w(:,ind.b),1); difsubv=mean(subhs.s(:,ind.b).^2,1);
difrealw=mean(realhs.w(:,ind.b),1); difrealv=mean(realhs.s(:,ind.b).^2,1);
    
% calculate the factors between real and sub weights and variances
% for pair sum
fws=sumrealw./sumsubw;
fvs=sumrealv./sumsubv;
% and for pair diff
fwd=difrealw./difsubw;
fvd=difrealv./difsubv;


for sdir=1:size(subac,2)  % over scan direction

  for j=1:size(subac,1) % over detector pairs

    %  RE-WEIGHT PAIR SUM QUANTITIES
    
    ac.wsum=subac(j,sdir).wsum.*fws(j);
    ac.wz=subac(j,sdir).wz.*fws(j);
    ac.wwv=subac(j,sdir).wwv.*fws(j).^2.*fvs(j);
    ac.switime=subac(j,sdir).switime.*fws(j);
    
    % RE-WEIGHT PAIR DIFF QUANTITIES
    
    ac.w=subac(j,sdir).w.*fwd(j);
    ac.wcz=subac(j,sdir).wcz*fwd(j);
    ac.wsz=subac(j,sdir).wsz*fwd(j);
    ac.wcc=subac(j,sdir).wcc*fwd(j);
    ac.wss=subac(j,sdir).wss*fwd(j);
    ac.wcs=subac(j,sdir).wcs*fwd(j);
    ac.wwccv=subac(j,sdir).wwccv.*fwd(j).^2.*fvd(j);
    ac.wwssv=subac(j,sdir).wwssv.*fwd(j).^2.*fvd(j);
    ac.wwcsv=subac(j,sdir).wwcsv.*fwd(j).^2.*fvd(j);
    % wzdiff is only needed when making coaddtype=3 pair diff maps
    % but is accumulated in reduc_makepairmaps as of 07/02/13
    ac.wzdiff=subac(j,sdir).wzdiff.*fwd(j);
    ac.dwitime=subac(j,sdir).dwitime.*fwd(j);

    if(isfield(subac(j,sdir),'wcd')) 
      for k=1:size(subac(j,sdir).wcd,2)
        ac.wcd{k}=subac(j,sdir).wcd{k}.*fwd(j);
        ac.wsd{k}=subac(j,sdir).wsd{k}.*fwd(j);
      end
    end

  end % over detector pairs
end % over scan direction

% these fields do not include weights, or they are removed from simulations in
% reduc_coaddpairmas anyway, so can keep them from the flavor substitution ac structure 
[ac.sitime]=subac.sitime;
[ac.swmax]=subac.swmax;
[ac.ditime]=subac.ditime;
[ac.dwmax]=subac.dwmax;

% reorder the fields so they are in the same order as typical
ac=orderfields(ac,subac);

return
