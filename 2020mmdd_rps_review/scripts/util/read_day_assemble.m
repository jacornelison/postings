function d=read_day_assemble(tag,acdc)
% d=read_day_assemble(tag,acdc)
%
% assemble "raw" files to make complete data structure
% (Matlab cant handle .mat files > 2GB but can hold bigger
%  structures in memory if available.)
%
% acdc='dc' means add back in DC level (staying in units of high gain volts)

if(~exist('acdc','var'))
  acdc='ac';
end

load(sprintf('data/raw_%s',tag));

load(sprintf('data/raw_%s_adc%d',tag,1));
d.lockin.adcData=x; clear x

load(sprintf('data/raw_%s_adc%d',tag,2));
d.lockin.adcData=[d.lockin.adcData,x]; clear x

load(sprintf('data/raw_%s_adc%d',tag,3));
d.lockin.adcData=[d.lockin.adcData,x]; clear x

load(sprintf('data/raw_%s_adc%d',tag,4));
d.lockin.adcData=[d.lockin.adcData,x]; clear x

% add in DC level if requested
if(strcmp(acdc,'dc'))
  dc=single(d.lockin.dcLevel);
  % add dc offset only when in high gain mode
  dc(d.lockin.outputMode~=1,:)=0;
  dc=repmat(shiftdim(dc,-1),[100,1,1]);
  dc=reshape(dc,[100*size(d.lockin.dcLevel,1),size(d.lockin.dcLevel,2)]);
  % add in dc offset in high gain units
  % see email from Jamie 10 Aug 2005 for explaination of scale factor
  d.lockin.adcData=d.lockin.adcData+200*(4.932/4096)*dc;
  clear dc
  % when in low gain mode scale to high gain units
  mode=repmat((d.lockin.outputMode~=1)',[100,1]);
  mode=cvec(mode);
  mode=repmat(mode,[1,96]);
  d.lockin.adcData(mode)=d.lockin.adcData(mode)*100;
  clear mode
end

return

d.t=d.frame.utc(:,2)/1e3;
d.tf=d.lockin.mjdMs/1e3;
d.td=d.frame.utc(:,1)+d.frame.utc(:,2)/86400e3;

subplot(3,1,1)
plot(d.t,d.frame.features)
subplot(3,1,2)
plot(d.tf,d.lockin.adcData(:,2))
subplot(3,1,3)
plot(d.t,d.lockin.dcLevel(:,2))
