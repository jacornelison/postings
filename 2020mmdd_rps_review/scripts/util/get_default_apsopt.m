function apsopt=get_default_apsopt(apsopt)
% apsopt=get_default_apsopt(apsopt)
if(~exist('apsopt','var'))
  apsopt=[];
end

if(~isfield(apsopt,'pure_b'))
  apsopt.pure_b='normal';
end
if(~isfield(apsopt,'update'))
  apsopt.update=0;
end
if(~isfield(apsopt,'howtojack'))
  apsopt.howtojack='dim2';
end
if(~isfield(apsopt,'ukpervolt'))
  apsopt.ukpervolt={get_ukpervolt()};
end
if(~isfield(apsopt,'overrideukpervolt'))
  apsopt.overrideukpervolt=0;
end
if(~iscell(apsopt.ukpervolt))
  if length(apsopt.ukpervolt==1)
    apsopt.ukpervolt={apsopt.ukpervolt};
  elseif length(apsopt.ukpervolt==2)
    apsopt.ukpervolt={apsopt.ukpervolt(1),apsopt.ukpervolt(2)};
  end
end
if(~isfield(apsopt,'hostdir'))
  apsopt.hostdir = {''};
end
if(~isfield(apsopt,'commonmask'))
  apsopt.commonmask = false;
end
if(~isfield(apsopt,'commonmasksel'))
  apsopt.commonmasksel = [];
end
if(~isfield(apsopt,'coaddrx'))
  apsopt.coaddrx = false;
end
if(~isfield(apsopt,'save_coaddopts'))
  apsopt.save_coaddopts = true;
end
if(~isfield(apsopt,'bintype'))
  apsopt.bintype = 'bicep_norm';
end
if(~isfield(apsopt,'makebpwf'))
  apsopt.makebpwf = 0;
end
if(~isfield(apsopt,'overall'))
  apsopt.overall = false;
end
if(~isfield(apsopt,'intorx'))
  apsopt.intorx = false;
end
if(~isfield(apsopt,'polrot'))
  apsopt.polrot = [0];
end
if(~isfield(apsopt,'purifmatname'))
  apsopt.purifmatname = {};
end
if(~isfield(apsopt,'pure_e'))
  apsopt.pure_e = 0;
end
if(~isfield(apsopt,'smoothvarmaps'))
  apsopt.smoothvarmaps=1;
end
if(~isfield(apsopt,'random_order'))
  apsopt.random_order=0;
end
if(~isfield(apsopt,'pad_map'))
  apsopt.pad_map=true;
end

return
