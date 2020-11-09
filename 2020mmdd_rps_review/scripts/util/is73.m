function ind=is73(files)
% is73(files)
%
% For input cell array files, returns logical array that is true of file is a version
% 7.3 matfile and false otherwise.
%
% ex. ind=is73({'a.mat','b.mat'});

if ~iscell(files)
  files={files};
end

ind=false(size(files));

for i=1:numel(files)
  [dum,f]=system(sprintf('file %s',files{i}));
  x=strfind(f,'mat-file');
  if ~isempty(x)
    x=strfind(f,'0x0200');
    if ~isempty(x)
      ind(i)=true;
    end
  end
end

return