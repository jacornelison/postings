function d=downsamp_fastreg(d,n,lf)
% d=downsamp_fastreg(d,n,lf)
%
% Downsample all the fast registers
%
% Modified to meet apparent phasing of fast/slow regs in BICEP data by
% making time.utcfast and time.utcslow regs match

disp('downsamp_fastreg...')

% search through d structure for fast sample registers and downsample
% them

names=fieldnames(d);

for i=1:length(names)
  if isstruct(d.(names{i}))
    d.(names{i})=downsamp_fastreg(d.(names{i}),n,lf);
  else
    if size(d.(names{i}),1)==lf
      d.(names{i})=downsample(d.(names{i}),n);
    end
  end
end

return
