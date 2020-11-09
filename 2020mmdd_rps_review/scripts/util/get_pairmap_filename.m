function fname=get_pairmap_filename(tags,mapopt)
% fname=get_pairmap_filename(tag,mapopt)
%
% Return the path+filename for a pairmap as defined by the tag and mapopt.
% 
% tags is either a string or a cell array of strings. In the latter case, the output is
% a cell array of strings.
%
% mapopt must have the following fields defined:
%  - sernum
%  - weight
%  - gs
%  - proj
%  - filt
%  - deproj
%
% ex. mapopt=get_default_mapopt;
%     fname=get_pairmap_filename(mapopt);

% If input is not a cell array, make it so
if ~iscell(tags)
  tags={tags};
  wascell=0;
else
  wascell=1;
end

for i=1:numel(tags)

  tag=tags{i};
  
  % append _gs if ground sub on
  if(mapopt.gs==1)
    gs='_gs';
  else
    gs='';
  end
  
  % append proj type if non-standard
  if(~strcmp(mapopt.proj,'radec'))
    proj=['_',mapopt.proj];
  else
    proj='';
  end
  
  % concat sum/diff filts for output filename
  if(iscell(mapopt.filt))
    filt=[mapopt.filt{:}];
  else
    filt=mapopt.filt;
  end
  
  % append deprojection type
  if iscell(mapopt.deproj)
    deproj = '_dp';
  else
    if numel(mapopt.deproj)>1 || any(mapopt.deproj)
      deproj='_dp';
    else
      deproj='';
    end
  end
  
  % make the filename
  fname{i}=sprintf('pairmaps/%s/%s/%s_filt%s_weight%1d%s%s%s.mat',...
                mapopt.sernum(1:4),...
                mapopt.sernum(5:end),...
                tag,...
                filt,...
                mapopt.weight,...
                gs,...
                proj,...
                deproj);
end

if ~wascell
  fname=fname{1};
end

return

