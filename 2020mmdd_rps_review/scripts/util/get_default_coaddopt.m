function coaddopt=get_default_coaddopt(coaddopt)
% coaddopt=get_default_coaddopt(coaddopt)
%
% Fills in the coaddopt structure with defaults for use with reduc_coaddpairmaps.
%
% If coaddopt.sernum is not specified on input, coaddopt.sernum and coaddopt.realpairmapset
% are not defined on output.
%
% Coaddopt fields common to mapopt (as defined in get_default_mapopt) really should
% return the same values for defaults, although they don't have to.
% 
% ex. coaddopt.filt='p0';
%     coaddopt.sernum='12050011';
%     coaddopt=get_default_coaddopt(coaddopt);
%               OR
%     coaddopt=get_default_coaddopt;


if(~exist('coaddopt','var'))
  coaddopt=[];
end

if(~isfield(coaddopt,'coaddtype'))
  coaddopt.coaddtype=0;
end
if(~isfield(coaddopt,'jacktype'))
  switch coaddopt.coaddtype
    case {0,1}
      coaddopt.jacktype='0123456789abcde'; % make all jackknives...
    otherwise
      coaddopt.jacktype='0'; % ...unless coaddtype is per channel
  end
elseif ~isstr(coaddopt.jacktype)
  % Convert numeric vector to string
  for k=1:numel(coaddopt.jacktype)
    x(k)=num2str(coaddopt.jacktype(k));
  end
  coaddopt.jacktype=x;
end
  
if(~isfield(coaddopt,'filt'))
  coaddopt.filt='p3';
end
if(~isfield(coaddopt,'gs'))
  coaddopt.gs=1; % default is now gs on
end
if(~isfield(coaddopt,'weight'))
  coaddopt.weight=3; % r-var-post-scanset
end
if(~isfield(coaddopt,'chflags'))
  coaddopt.chflags=get_default_chflags;
end
if(~isfield(coaddopt,'proj'))
  coaddopt.proj='radec';
end
if(~isfield(coaddopt,'sernum'))
  disp(['No serial number (coaddopt.sernum) specified. '...
           'coaddopt.sernum and coaddopt.realpairmapset will not be defined.'])
elseif ~isfield(coaddopt,'realpairmapset')
  % this is an "evil" option that could lead to confusion - don't
  % use it!
  coaddopt.realpairmapset = coaddopt.sernum(1:4);
end
if(~isfield(coaddopt,'cut'))
  % if not specified get standard cuts
  %if given tags, use the first one for the year
  if(isfield(coaddopt,'tags'))
    coaddopt.cut=get_default_round2_cuts([],coaddopt.tags{1}(1:4));
  else
    coaddopt.cut=get_default_round2_cuts;
  end
end
if(~isfield(coaddopt,'daughter'))
  coaddopt.daughter='a';
end
if(~isfield(coaddopt,'deproj'))
  coaddopt.deproj=false;
end
if(~isfield(coaddopt,'deproj_timescale'))
  coaddopt.deproj_timescale='by_phase';
end
if(~isfield(coaddopt,'deproj_src'))
  coaddopt.deproj_src=[];
end
if(~isfield(coaddopt,'save_cuts'))
  coaddopt.save_cuts=true;
end
if(~isfield(coaddopt,'do_covariance'))
  coaddopt.do_covariance=0;
end
if(~isfield(coaddopt,'sign_flip_type'))
  coaddopt.sign_flip_type=false;
end

return
