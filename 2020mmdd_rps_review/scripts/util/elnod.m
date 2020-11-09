function [eng,dat]=elnod(d,p,en,ch)
% [eng,dat]=elnod(d,p,en,ch)
%
%   Extract rel cal factors from elnods.
% 
%   Outputs:
%     eng contains calculated relative gains:
%       [offset,gain,goodness-of-fit,nan-count];
%     dat contains el nod data;

disp('elnod...')

% ensure full size output regardless of which channels actually get
% fit - 2 fit pars plus gof
n=length(en.s); m=size(d.mce0.data.fb,2);

eng=NaN*zeros(n,m,4);
dat(n,m).x=[]; dat(n,m).y=[]; dat(n,m).m=[];

for i=1:length(en.s)
  
  % fetch the data
  s=en.sf(i); e=en.ef(i);
  npts=length(s:e);
  y=double(d.mce0.data.fb(s:e,:));
  az=d.pointing.hor.az(s:e);
  el=d.pointing.hor.el(s:e);
  dk=d.pointing.hor.dk(s:e);
  
  for j=ch
    
    % calculate az/el for this feed
    % (By comparison to "known working" line in reduc_makecomap I
    % believe this is right - since we are working with elevation
    % rather than Dec we want to offset in the opposite direction
    % meaning 180 deg rotation of bearing angle.) 
    % Note this is in principle wrong: az/el and ra/dec are
    % opposite handed coordinate systems, so we should be using
    % -(p.theta-90-dk).  This will result in a mirror flip,
    % which should only affect the resulting "dummy" az, so
    % I'm leaving it for now.  -- RWO 2013-apr-30
    [elf,dummy]=reckon(el,az,p.r(j),p.theta(j)-90-dk);

    za=90-elf;
    
    am=1./cos(za*pi/180);
    X=[ones(size(am)) am];

    if sum(~isnan(y(:,j)))>=10
      b=regress(y(:,j),X);
    
      % calc chisq as some estimate of the goodness of fit
      
      gof=nansum((y(:,j)-X*b).^2)/npts/nanstd(diff(y(:,j))).^2;;
    else
      % All or most samples are NaN - fill output with NaN
      b=[NaN; NaN];
      gof=NaN;
    end

    % calc number of NaN / dropped / deglitched samples
    % as a parameter for quality cuts
    nnan=sum(isnan(y(:,j)));

    % store the fit params and the gof
    eng(i,j,:)=[b;gof;nnan];
      
    % store the data which was used in the fit
    dat(i,j).x=am; dat(i,j).y=y(:,j); dat(i,j).m=X*b;
    dat(i,j).el=elf;

    if(0)
      clf
      plot(am,y(:,j));
      hold on; plot(am,polyval2(b,am),'r'); hold off
      xlabel('air masses'); ylabel('ADC (V)');
      title(sprintf('fitting %s',p.channel_name{j}));
      axis tight
      pause
    end
  end
end

return
