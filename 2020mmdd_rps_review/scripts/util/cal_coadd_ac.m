function ac=cal_coadd_ac(ac,ukpervolt,coaddopt,scaleV,scaleS)
% ac=cal_coadd_ac(ac,ukpervolt,coaddopt,NV,NS)
%
% Apply scale factor (ukpervolt) to coadd acs
%
% coaddopt is needed to pick the proper weighting scheme
% if not handed in, the regular weight=3 is assumed
%
% scaleV:  (1 default) scale the variance portions of the ac, or not (0)
% scaleS:  (1 default) scale the signal portions of the ac, or not (0)
%
% for regular application of ukpervolt both signal 
% and variances are scaled scaleV=scaleS=1 (default)
%
% for fixing abscal issues in the sims, or applying a scaling
% factor such as r=0.1 -> r=0.2 one wants to scale the signal
% but leave the variances alone: scaleV=0, scaleS=1
%  
% Note: For matrix the selection of the portion is not implemented
% it will realize the default case scaleV=scaleS=1

% Assume weight=3 if not specified
if ~exist('coaddopt','var')||isempty(coaddopt)
  coaddopt.weight=3;
end

if ~exist('scaleV','var')||isempty(scaleV)
  scaleV=1;
end

if ~exist('scaleS','var')||isempty(scaleS)
  scaleS=1;
end

% coaddopt from maps combined over different observations
if iscell(coaddopt)
  coaddopt=coaddopt{1};
else
  coaddopt=coaddopt;
end

% select the exponent multipliers for the scaling explicitly:
% 1. the weighting scheme
if coaddopt.weight==0
  NW=0; % uniform
elseif any(coaddopt.weight==[1,2,3])
  NW=1; % inverse variance
end
% 2. the variance scaling:
NV=1;
if ~scaleV
  NV=0;
  % if the variance is not scaled also don't scale the weights
  NW=0;
end
% 3. the signal scaling:
NS=1;
if ~scaleS
  NS=0;
end

%combcomap just concatenates ac, so have to loop through those
if iscell(ac)
  for kk=1:length(ac)
    ac{kk}=cal_ac(ac{kk},ukpervolt,NV,NS,NW);
  end
else
  ac=cal_ac(ac,ukpervolt,NV,NS,NW);
end

return


%%%%%%%%%%%%%%%%%%%%%%%
% breakout function
function ac=cal_ac(ac,ukpervolt,NV,NS,NW)

% If ukpervolt is single valued, expand to size of ac first dimension 
if(size(ac,1)>1 && numel(ukpervolt)==1)
  ukpervolt=repmat(ukpervolt,[size(ac,1),1]);
end

% apply the cal factors
for i=1:size(ac,1)
  for j=1:size(ac,2)
    % pair sum
    if isfield(ac(i,j),'wsum')
      ac(i,j).wsum=ac(i,j).wsum*ukpervolt(i)^(-2*NW);
    end
    if isfield(ac(i,j),'wz')
      ac(i,j).wz=ac(i,j).wz*ukpervolt(i)^(-2*NW+NS);
    end
    if isfield(ac(i,j),'wwv')
      ac(i,j).wwv=ac(i,j).wwv*ukpervolt(i)^(-4*NW+2*NV);
    end
    if isfield(ac(i,j),'switime')
      ac(i,j).switime=ac(i,j).switime*ukpervolt(i)^(-2*NW);
    end
    if isfield(ac(i,j),'swmax')
      ac(i,j).swmax=ac(i,j).swmax*ukpervolt(i)^(-2*NW);
    end
    % pair diff
    if isfield(ac(i,j),'w')
      ac(i,j).w=ac(i,j).w*ukpervolt(i)^(-2*NW);
    end
    if isfield(ac(i,j),'wcz')
      ac(i,j).wcz=ac(i,j).wcz*ukpervolt(i)^(-2*NW+NS);
      ac(i,j).wsz=ac(i,j).wsz*ukpervolt(i)^(-2*NW+NS);
      ac(i,j).wcc=ac(i,j).wcc*ukpervolt(i)^(-2*NW);
      ac(i,j).wss=ac(i,j).wss*ukpervolt(i)^(-2*NW);
      ac(i,j).wcs=ac(i,j).wcs*ukpervolt(i)^(-2*NW);
    end
    if isfield(ac(i,j),'wwccv')
      ac(i,j).wwccv=ac(i,j).wwccv*ukpervolt(i)^(-4*NW+2*NV);
      ac(i,j).wwssv=ac(i,j).wwssv*ukpervolt(i)^(-4*NW+2*NV);
      ac(i,j).wwcsv=ac(i,j).wwcsv*ukpervolt(i)^(-4*NW+2*NV);
    end
    if isfield(ac(i,j),'wzdiff')
      ac(i,j).wzdiff=ac(i,j).wzdiff*ukpervolt(i)^(-2*NW+NS);
    end
    if isfield(ac(i,j),'dwitime')
      ac(i,j).dwitime=ac(i,j).dwitime*ukpervolt(i)^(-2*NW);
    end
    if isfield(ac(i,j),'dwmax')
      ac(i,j).dwmax=ac(i,j).dwmax*ukpervolt(i)^(-2*NW);
    end
    % poly-sub and gnd-sub maps
    if isfield(ac(i,j),'wz_psub')
      ac(i,j).wz_psub=ac(i,j).wz_psub*ukpervolt(i)^(-2*NW+NS);
      ac(i,j).wz_gsub=ac(i,j).wz_gsub*ukpervolt(i)^(-2*NW+NS);
    end
    if isfield(ac(i,j),'wcz_psub')
      ac(i,j).wcz_psub=ac(i,j).wcz_psub*ukpervolt(i)^(-2*NW+NS);
      ac(i,j).wsz_psub=ac(i,j).wsz_psub*ukpervolt(i)^(-2*NW+NS);
      ac(i,j).wcz_gsub=ac(i,j).wcz_gsub*ukpervolt(i)^(-2*NW+NS);
      ac(i,j).wsz_gsub=ac(i,j).wsz_gsub*ukpervolt(i)^(-2*NW+NS);
    end
    % deprojection templates
    if isfield(ac(i,j),'wcd')
      for k=1:numel(ac(i,j).wcd)
        ac(i,j).wcd{k}=ac(i,j).wcd{k}*ukpervolt(i)^(-2*NW+NS);
        ac(i,j).wsd{k}=ac(i,j).wsd{k}*ukpervolt(i)^(-2*NW+NS);
      end
    end

    %matrix pieces
    if isfield(ac(i,j),'awa')
      ac(i,j).awa=ac(i,j).awa*ukpervolt(i)^(-2*NW);
    end
    if isfield(ac(i,j),'awwva')
      ac(i,j).awwva=ac(i,j).awwva*ukpervolt(i)^(-4*NW+2*NV);
    end
    if isfield(ac(i,j),'awta')
      ac(i,j).awta=ac(i,j).awta*ukpervolt(i)^(-2*NW);
    end
    if isfield(ac(i,j),'awssa')
      ac(i,j).awssa=ac(i,j).awssa*ukpervolt(i)^(-2*NW);
      ac(i,j).awcca=ac(i,j).awcca*ukpervolt(i)^(-2*NW);
      ac(i,j).awcsa=ac(i,j).awcsa*ukpervolt(i)^(-2*NW);
    end
    if isfield(ac(i,j),'awdiffa')
      ac(i,j).awdiffa=ac(i,j).awdiffa*ukpervolt(i)^(-2*NW);
    end
    if isfield(ac(i,j),'awwssva')
      ac(i,j).awwssva=ac(i,j).awwssva*ukpervolt(i)^(-4*NW+2*NV);
      ac(i,j).awwccva=ac(i,j).awwccva*ukpervolt(i)^(-4*NW+2*NV);
      ac(i,j).awwcsva=ac(i,j).awwcsva*ukpervolt(i)^(-4*NW+2*NV);
    end
    if isfield(ac(i,j),'awpa')
      ac(i,j).awpa=ac(i,j).awpa*ukpervolt(i)^(-2*NW);
    end    
    if isfield(ac(i,j),'awgfha')
      ac(i,j).awgfha=ac(i,j).awgfha*ukpervolt(i)^(-2*NW);
      ac(i,j).awcgfcha_awcgfsha=ac(i,j).awcgfcha_awcgfsha*ukpervolt(i)^(-2*NW);
      ac(i,j).awsgfcha_awsgfsha=ac(i,j).awsgfcha_awsgfsha*ukpervolt(i)^(-2*NW);
    end

  end
end

return
