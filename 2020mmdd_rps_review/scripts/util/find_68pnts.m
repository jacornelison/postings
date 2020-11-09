function [cont68,lo68,hi68]=find_68pnts(pdf,vals,method,nsigma)
% function [cont68,lo68,hi68]=find_68pnts(pdf,vals,method)
% find contour of equal likelihood that encloses 68% of the pdf
% starting from the most likely point of the pdf
%  
% method    =   'interp': performs search using interpolation
%           =   'numeric': equivalent, performs search using
%               numeric integration, this one will survive double peaks
%               in the pdf. In the case of two peak the outermost
%               values will be chosen.
%
% nsigma    = allows to go away from 68%=+-1sigma to nsimga, default 1

if ~exist('nsigma','var')
  nsigma = 1;
end


if ~exist('method','var')
  method = 'interp'
end

[dum indx]=max(pdf);
maxlike=vals(indx);
sigma=0.3413*nsigma;

switch method
 case 'interp'
  normint=abs(cumsum(pdf)./sum(pdf));
  normlik=abs(pdf./sum(pdf));
  poslik=normlik(indx:end);
  neglik=normlik(indx:-1:1);
  rpos=vals(indx:end);
  rneg=vals(indx:-1:1);
  lowlim=1e-3;
  contvals=max(normlik):-.001:lowlim;
  for ii=1:length(contvals)
    n1=interp1(neglik(neglik>lowlim),rneg(neglik>lowlim),contvals(ii));
    p1=interp1(poslik(poslik>lowlim),rpos(poslik>lowlim),contvals(ii));
    p1x=interp1(vals,normint,p1);
    n1x=interp1(vals,normint,n1);   
    totint(ii)=p1x-n1x;
  end
  cont68=interp1(totint(~isnan(totint)),contvals(~isnan(totint)),sigma*2);
  lo68=interp1(neglik(neglik>lowlim),rneg(neglik>lowlim),cont68);
  hi68=interp1(poslik(poslik>lowlim),rpos(poslik>lowlim),cont68);
  
 case 'numeric'
  vals_f = linspace(vals(1),vals(end),100000);
  % make sure that the maxlike r value is part of the new range:
  %If there are multiple maxima, arbitrarily pick the first
  indx_f = find(abs(maxlike-vals_f)==min(abs(maxlike-vals_f)), 1, 'first');
  vals_f = vals_f+maxlike-vals_f(indx_f);
  % construct a finely binned pdf:
  normlik_f=interp1(vals,pdf,vals_f);
  normlik_f(isnan(normlik_f))=0;
  norm = sum(normlik_f);
  normlik_f=abs(normlik_f/norm);
  % now we integrate the ll function numerically in a range enclosed
  % by constant likehood. We start from the maximum likehood
  contvals = linspace(normlik_f(indx_f),0,1000);
  totint = zeros(size(contvals));
  rhi  = zeros(size(contvals));
  rlo  = zeros(size(contvals));
  for ii=1:length(contvals)
    binint = normlik_f>=contvals(ii);
    rhi(ii)= max(vals_f(binint));
    rlo(ii)= min(vals_f(binint));
    totint(ii) = sum(normlik_f(binint));
  end
  % the rebinning changes the absolute values of the pdf,
  % undo the nomalization and apply the normalization of
  % regular pdf:
  [ti,contvals]=pre_interp1(totint',contvals');
  [ti,rlo]=pre_interp1(totint',rlo');
  [ti,rhi]=pre_interp1(totint',rhi');
  cont68= interp1(ti,contvals,sigma*2)*norm/sum(pdf);
  lo68  = interp1(ti,rlo,sigma*2);
  hi68  = interp1(ti,rhi,sigma*2);

end

return

function [x_out,y_out]=pre_interp1(x_in,y_in)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 20130204 GPT
%
% The function interp1 often complains,
% "The grid vectors are not strictly monotonic increasing."
% This simple function sorts and rebins your x and y input and turns it
% into an output that Matlab can handle.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs:
%
%  x_in = column vector where data is sampled
%
%  y_in = column vector of data values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Outputs:
%
%  x_out = sorted, unique column vector where data is sampled
%
%  y_out = column vector of averaged data values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(x_in,2)>1 || size(y_in,2)>1
    x_out='please use a column vector for now';
    y_out=[];
    return
end

% ignore NaN
nan_index=union(find(isnan(x_in)),find(isnan(y_in)));
x_in(nan_index)=[];
y_in(nan_index)=[];

% Sort into ascending x.  Note that we can't use sortrows because the y
% data may be complex, which will confuse sorting real x.
[x_in,ix]=sort(x_in);
y_in=y_in(ix);

% If there are multiple y values for a given x, take their mean.
x_out=unique(x_in);
N=length(x_out);
y_out=nan(N,1);
n=hist(x_in,x_out);
n_sum=cumsum(n);
y_out(1)=mean(y_in(1:n_sum(1)));
for jj=2:N
    y_out(jj)=mean(y_in((n_sum(jj-1)+1):n_sum(jj)));
end

return
