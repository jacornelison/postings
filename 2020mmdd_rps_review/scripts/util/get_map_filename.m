function fname=get_map_filename(coaddopt)
% fname=get_map_filename(coaddopt)
%
% Return the path+filename for a coadded map as defined by mapopt and coaddopt.
% 
% coaddopt must have the following fields defined:
%  - sernum
%  - weight
%  - gs
%  - proj
%  - filt
%  - deproj
%  - daughter
%  - jacktype
%  - coaddtype
%
% If coaddopt.jacktype is an array, the list of filenames returned is a cell array

% append _gs if ground sub on
if(coaddopt.gs==1)
  gs='_gs';
else
  gs='';
end

% append proj type if non-standard
if(~strcmp(coaddopt.proj,'radec'))
  proj=['_',coaddopt.proj];
else
  proj='';
end

% concat sum/diff filts for output filename
if(iscell(coaddopt.filt))
  filt=[coaddopt.filt{:}];
else
  filt=coaddopt.filt;
end

% compose filename extension
filenameext=sprintf('filt%s_weight%1d%s%s',filt,coaddopt.weight,gs,proj);

% only append the second digit to jacktype if non default coaddtype
if(coaddopt.coaddtype>0)
  coaddtype=sprintf('%1d',coaddopt.coaddtype);
else
  coaddtype='';
end

for ii=1:numel(coaddopt.jacktype)  
  for j=1:size(coaddopt.deproj,1)

    % append deprojection type
    if numel(coaddopt.deproj(j,:))>1 || coaddopt.deproj(j,1)
      deproj='_dp';
      for i=1:size(coaddopt.deproj,2)
        deproj=[deproj,num2str(coaddopt.deproj(j,i))];
      end
    else
      deproj='';
    end

    % make the filename
    fname{ii,j}=sprintf('maps/%s/%s_%s_%s%s_jack%s%s.mat',...
                   coaddopt.sernum(1:4),...
                   coaddopt.sernum(5:end),...
                   coaddopt.daughter,...
                   filenameext, deproj,...
                   coaddopt.jacktype(ii),...
                   coaddtype);
    if isfield(coaddopt,'sign_flip_rlz') 
      % we are generating noise realizations from real data
      % give this the identifier 6 and change from real to rlz number
      fname{ii,j} = strrep(fname{ii,j},'real',sprintf('%03d6',coaddopt.sign_flip_rlz));
      % if applying to sim data make sure we at least capture the right rlz:
      fname{ii,j}(11:13)=sprintf('%03d',coaddopt.sign_flip_rlz);
    end
  end
end

if numel(fname)==1
  fname=fname{1};
end

return

