function r=calc_devs(r,chibins,meansim)
% r=calc_devs(r,chibins)
%
% calculate bandpower deviations given:
% r.expv - expectation values (often mean of sig+n sims)
% r.real - real bandpowers
% r.derr - bandpower errorbars (usually form mean of sig+n sims)
% r.sim - sim bandpowers
%
% Outputs:
% r.rdevs - real data bandpower deviations in units of sigma (from r.derr)
% r.sdevs - sim bandpower deviations in units of sigma
% r.rsumdev - sum devs over chibins
% r.ssumdev - sum devs over chibins
% r.ptet_sumdev - pte of rsumdev calc against theory distribution
% r.ptes_sumdev - pte of rsumdev calc against sim distribution

if(~exist('meansim','var'))
  meansim=[];
end

if(meansim)
  % set expectation value to mean of s+n sims
  for j=1:length(r)
    r(j).expv=mean(r(j).sim,3);
  end
end

if(~isfield(r,'derr'))
  r=get_bpcov(r);
end

for j=1:length(r)
  for i=1:size(r(j).real,2);
   
    % model
    m=r(j).expv(:,i);
    
    % calc normalized deviation for real
    d=r(j).real(:,i);
    e=r(j).derr(:,i);
    r(j).rdevs(:,i)=(d-m)./e;
    
    % calc normalized deviation for sig+noi
    for k=1:size(r(j).sim,3)
      d=r(j).sim(:,i,k);
      r(j).sdevs(:,i,k)=(d-m)./e;
    end
    
    % calc sum of devs
    r(j).rsumdev(i)=sum(r(j).rdevs(chibins{i},i));
    r(j).ssumdev(i,:)=sum(r(j).sdevs(chibins{i},i,:),1);
  
    % calc pte sim
    r(j).ptes_sumdev(i)=sum(r(j).rsumdev(i)<r(j).ssumdev(i,:))/size(r(j).ssumdev,2);

    % Ignoring bandpower covariance the sumdevs is simply a sum of
    % length(chibins{i}) Gaussian random numbers each with unit std -
    % i.e. it should have std of 1/sqrt(length(chibins{i}))
    % However in the presence of covariance it will tend to be larger
    % - choose to use the fitted width and mean to allow for such
    % effect. This "theory" number is therefore basically an
    % extrapolation of the sim histogram assuming Gaussianity
    r(j).ptet_sumdev(i)=1-normcdf(r(j).rsumdev(i),mean(r(j).ssumdev(i,:)),std(r(j).ssumdev(i,:)));
  
    % copy chibins
    r(j).chibins{i}=chibins{i};
  
  end
end

return
