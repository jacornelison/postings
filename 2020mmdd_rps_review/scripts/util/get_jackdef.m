function [jknames freqjknum]=get_jackdef(jkret)
%
% This function should be used to control
% the naming conventions of the jackknives
%
% input:  optional jkret, return only jackname requested
% output: if no input, return all jackknife names
%         optional, return frequency jack#

freqjknum=[];
jknames={'Signal','Deck jackknife','Scan Dir jackknife','Tag Split jackknife',...
         'Tile jackknife','Phase jackknife','Mux Col jackknife',...
         'Alt Deck jackknife','Mux Row jackknife','Tile/Deck jackknife',...
         'Focal Plane inner/outer jackknife','Tile top/bottom jackknife',...
         'Tile inner/outer jackknife','Moon jackknife','A/B offset best/worst',...
         'Map jackknife','Spectral jackknife','Pulsetube sensitive most/least',...
         'Receiver jackknife'};
jkletters={'0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f','s','g','r'};
jacks = containers.Map(jkletters,jknames);

if ~exist('jkret') || isempty(jkret)
  jkret=jkletters(2:end-1);
end

% Convert to cell
if isstr(jkret)
  x=jkret;
  clear jkret
  for k=1:numel(x)
    jkret{k}=x(k);
  end
end

if ~iscell(jkret)
  % Must be a scalar or array of numbers
  for k=1:numel(jkret)
    y{k}=num2str(jkret(k));
  end
  jkret=y;
end

% Make sure everything is a string
for k=1:numel(jkret)
  jkret{k}=num2str(jkret{k});
end

% Output the jack names
for ii=1:length(jkret)
  k=strfind(char(jkletters)',jkret{ii});
  jknm{ii}=jknames{k};
end
jknames=jknm;

if numel(jknames)==1
  jknames=jknames{1};
end

return



