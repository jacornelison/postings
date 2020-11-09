function calib2csv(expt,yr,workdir,csvfname)

if isempty(expt)
  expt=get_experiment_name();
end

oldpwd=pwd;
cd(workdir)

% Year or tag date
dn=[];
if ischar(yr)
  if length(yr)>4
    dn=yr;
    yr=str2num(yr(1:4));
  else
    dn=[yr '0601'];
    yr=str2num(yr);
  end
else
  if yr>10000
    dn=num2str(yr);
    yr=floor(yr/10000);
  else
    dn=[num2str(yr) '0601'];
    yr=yr;
  end
end

% Read in the calibration files for the MCE's
if strcmp(expt,'bicep3')
  calibfunc = 'calib_bicep3_run8';
elseif strcmp(expt,'keck')
  if yr==2015
    calibfunc = 'calib_keck_pole2015';
  elseif yr==2013
    calibfunc = 'calib_keck_pole2013';
  elseif yr==2012
    calibfunc = 'calib_keck_pole2012';
  elseif yr==2011
    calibfunc = 'calib_keck_pole2011';
  end
else
  calibfunc = 'calib_bicep2_run8';
end
calib=eval(calibfunc);

% Output file name
if ~exist('csvfname','var') || isempty(csvfname)
  csvfname=fullfile('aux_data','detparams.csv');
end

% get_array_info so we know num channels, etc.
[p0,ind0]=get_array_info(dn);

% New CSV header
clear k
k.comments={['# Detector parameters for ' expt]};
if isempty(calibfunc)
  k.comments{2}=['# Calibration information passed as structure'];
else
  k.comments{2}=['# Calibration information from ' calibfunc];
end
k.comments{3}=['# Created ' getenv('USER') ', ' datestr(now)];
[tmpa tmpb tmpc]=fileparts(csvfname);
k.filename=[tmpb tmpc];
k.created=[datestr(now,'yyyymmdd') ' ' getenv('USER')];
k.fields={'r_sh','r_wire','r_bias','bits_bias','v_b_max','r_fb1','bits_fb1','v_fb1_max','m_fb1'};
k.units={'Ohm','Ohm','Ohm','(none)','V','Ohm','(none)','V','(none)'};
k.formats={'float','float','float','integer','float','float','integer','float','float'};

% Expand values to full chan list
clear p
p.gcp=p0.gcp;
for i=1:length(k.fields)
  cval = calib.(upper(k.fields{i}));
  p.(k.fields{i})=expand_to_detlist(cval,p0);
end

% And write out
ParameterWrite(csvfname,p,k);

cd(oldpwd)

return

function v=expand_to_detlist(v,p)

% Empty -- fill with zeros
if length(v)==0
  v=0.0*ones(size(p.gcp));
  return
end
% Single value -- expand
if length(v)==1
  v=v*ones(size(p.gcp));
  return
end
% Already full length -- keep as is
if length(v)==length(p.gcp)
  v=reshape(v,size(p.gcp));
  return
end
% One per MCE column as in BICEP3, same for every MCE
%  -- makes no sense, but expand as given
if length(v)==1+max(p.mce_col)
  tmp=zeros(size(p.gcp));
  for i=0:max(p.mce_col)
    tmp(p.mce_col==i)=v(i+1);
  end
  v=tmp;
  return
end
% One per unique MCE column as in Keck
if length(v)==(1+max(p.mce))*(1+max(p.mce_col))
  tmp=zeros(size(p.gcp));
  v=reshape(v,[],1+max(p.mce));
  for i=0:max(p.mce)
    for j=0:max(p.mce_col)
      tmp(p.mce==i & p.mce_col==j)=v(j+1,i+1);
    end
  end
  v=tmp;
  return
end
% Messed up wrong length but constant as for some
% fields in Keck -- just expand out constant value
% and don't worry about it
if length(unique(v)==1)
  v=unique(v)*ones(size(p.gcp));
  return
end

error(['Don''t know what to do with calib field of length ' num2str(length(v))]);

return

