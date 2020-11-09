function map=insert_common_mask(m,map,opt,mapsel)
% map=insert_common_mask(m,map,opt,mapsel)
%
% If opt='gmean' take the gmean mask from the maps specified
% in mapsel and copy to all. 
% If opt=integer copy the mask from that map to all.
%  
% When the product spectra are to be jointly considered (e.g. is a
% multispectral fit) then this is probably what we want. Note that
% it's not a total solution - because of filtering the modal
% coverage will still not be exactly the same.

% don't do this for jack maps - for some like tile the B2 and Keck
% coverage is radically different - don't want to eat this
% difference simply to make the poorly motivated jackknife cross
% spectrum better behaved
if(size(map,2)>1)
  % have this a warning instead of an error since the delta maps stacks
  % use this dimension.
  warning('commonmask for jack~=0 a funny thing to do');
end
% number of maps
n1=size(map,1);
na=numel(map);

% if no selection of maps is given to use for the gmean,
% use all in dimension 1:
if ~exist('mapsel','var') || isempty(mapsel)
  mapsel = 1:n1;
end

if ~isstr(opt)
  
  disp(sprintf('taking common mask from map %d',opt));
  Tw=map(opt).Tw;
  Pw=map(opt).Pw;
  
else
  switch opt
   case 'gmean'
    
    disp('inserting common mask as gmean of maps');

    % take product of masks
    Tw = map(mapsel(1)).Tw;
    Pw = map(mapsel(1)).Pw;
    for i=2:length(mapsel)
      Tw=Tw.*map(mapsel(i)).Tw;
      Pw=Pw.*map(mapsel(i)).Pw;
    end
  
    % take the n'th root
    Tw=Tw.^(1/length(mapsel));
    Pw=Pw.^(1/length(mapsel));
  end
end
  
% copy to all
for i=1:na
  map(i).Tw=Tw;
  map(i).Pw=Pw;
end

return
