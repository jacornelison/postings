function apsopt = aps_load_purifmat(apsopt)
% function apsopt = aps_load_purifmat(apsopt)
% load the purification matrix specified in apsopt.purifmatname
% into apsopt.purifmat. Avoids reloading the matrix for duplicate
% entries.

% if using matrix purification load the matrix once now
if ~isempty(apsopt.purifmatname)
  
  if ~iscell(apsopt.purifmatname)
    % if purifmatname is a simple string then assume the same one
    % is to be used for all maps
    purifmatname=apsopt.purifmatname;
    apsopt=rmfield(apsopt,'purifmatname');
    for i=1:length(apsopt.polrot)
      apsopt.purifmatname{i}=purifmatname;
    end
  end
  
  % fetch matrix for each entries in the cell:
  for ii=1:length(apsopt.purifmatname)

    % empty cell entries need to be passed forward:
    if isempty(apsopt.purifmatname{ii})
      apsopt.purifmat{ii}=[];
      continue;
    end

    % check if several maps get the same matrix, to avoid loading the file serval times:
    is_duplicate=0;
    for jj=1:ii-1
      if strcmp(apsopt.purifmatname{ii},apsopt.purifmatname{jj})
        disp(['Reuse purification matrix : ',apsopt.purifmatname{jj}])
        apsopt.purifmat{ii}=apsopt.purifmat{jj};
        is_duplicate=1;
        break
      end
    end

    if ~is_duplicate
      disp(['Loading purification matrix : ',apsopt.purifmatname{ii}])
      matfields = who('-file', apsopt.purifmatname{ii});
      if ~ismember('projmatopt', matfields)
        vars = {'m', 'obs_pixels', 'reob', 'pb'};
      else
        vars = {'m', 'projmatopt', 'pb'};
      end
      if apsopt.pure_e
        vars = [vars 'pe'];
      end
      apsopt.purifmat{ii} = load(apsopt.purifmatname{ii}, vars{:});
    end
  
  end % for

end

return
