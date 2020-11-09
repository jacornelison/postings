function fs=filter_fs(fs,filt)
% fs=filter_fs(fs,filt)
%
% filter scan pointers

% remove the windback scan from fs

ind=false(size(fs.s));
for i=1:length(fs.s)
  if(all(filt(fs.s(i):fs.e(i))))
    ind(i)=1;
  end
end

fs=structcut(fs,ind);

return
