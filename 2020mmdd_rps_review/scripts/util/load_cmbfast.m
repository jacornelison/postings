function mod=load_cmbfast(filename,coaddopt)
% mod=load_cmbfast(filename,coaddopt)
%
% Load cmbfast spectrum and normalize to uk^2
% Note Cs_l is "C script l"
%
% If arg coaddopt present will mod the spectra to reflect what was
% actually used in sim
  
if isempty(filename)
  filename='input_maps/official_cl/camb_wmap5yr_noB.fits';
end

disp(sprintf('reading cmb spectrum from %s',filename));

if strcmp(filename(end-3:end), 'fits')
  temp=(fitsread(filename,'table'));
  ncol=size(temp,2);
  mod.l=(0:length(temp{1})-1)';
  mod.Cs_l=zeros(length(temp{1}),6);
  for i=1:min([4,ncol])
    mod.Cs_l(:,i)=temp{i}.*mod.l.*(mod.l+1)*(1e6^2)/(2*pi);
  end
elseif strcmp(filename(end-3:end), '.txt')
  temp=load(filename);
  mod.l=temp(:,1);
  mod.Cs_l=(2.73^2)*(1e6^2)*[temp(:,2:end)];
else
  error('feed me a .fits or a .txt filename');
end
 
% if BB not in file set it to 0
if(size(mod.Cs_l,2)==3)
  mod.Cs_l(:,4)=0;
end

% re-order to make col TT,TE,EE,BB as Cynthia uses.
mod.Cs_l=mod.Cs_l(:,[1,4,2,3]);

% add TB,EB
mod.Cs_l(:,5:6)=0;

% add ET and BT, BE
mod.Cs_l(:,7)=mod.Cs_l(:,2);
mod.Cs_l(:,8:9)=0;


% scale from l(l+1)/2pi form to straight C_l
lcol=repmat(mod.l,[1,size(mod.Cs_l,2)]);
mod.C_l=(mod.Cs_l*2*pi)./(lcol.*(lcol+1));

% if arg coaddopt provided
if(exist('coaddopt','var'))
  
  % if was a sim run modify the spectra to reflect what really went in
  if(isfield(coaddopt.mapopt{1},'simopt'))
    
    % changed the name of this field - maintain backwards compat
    if(isfield(coaddopt.mapopt{1}.simopt,'sig'))
      sig=coaddopt.mapopt{1}.simopt.sig;
    elseif(isfield(coaddopt.mapopt{1}.simopt,'cmb'))
      sig=coaddopt.mapopt{1}.simopt.cmb;
    else
      %do nothing
    end
    
    switch sig
      case 'nopol'
	mod.Cs_l(:,2:4)=0;
      case 'onlypol'
	mod.Cs_l(:,1)=0;
      case 'none'
	mod.Cs_l(:,:)=0;
    end
  end
  
  % backwards compatibility catch - removed later
  if(~ischar(coaddopt.jacktype))
    coaddopt.jacktype=num2str(coaddopt.jacktype);
  end
  
  % if was a jackknife run zero the expected spectra
  if(coaddopt.jacktype~='0')
    mod.Cs_l(:,:)=0;
  end
  
end

return
