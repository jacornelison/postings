function ind=make_ind(p)
% ind=make_ind(p)
%
% generate convenient ind struture of index lists into array data
% structure for referencing bolometers of various classes
%

% Enumerate frequencies and receivers
freqs=unique(p.band(p.band~=0 & ~isnan(p.band)));
fstr=arrayfunc(@(f) sprintf('%03d',f), freqs);
rxs=unique(p.rx);

% Pre-determine the light, no-flag, etc masks since they're used many times
lmask=strcmp(p.type,'L');
amask=strcmp(p.pol,'A');

p.flag(isnan(p.flag))=0; % empty field reads in as NaN
fmask=(p.flag==0);


% channel type groups
ind.e=[1:length(p.band)]';      % everything
ind.l=find(lmask); % light
ind.d=find(strcmp(p.type,'D')); % dark
ind.o=find(strcmp(p.type,'O')); % open
for ff=1:length(freqs)
  ind.(['l' fstr{ff}])=find(lmask & p.band==freqs(ff));
end

% channels listed as good in pars file
ind.ge=find(fmask);  % demand no flag to count channel as good
ind.gl=find(lmask & fmask);
ind.gd=find(strcmp(p.type,'D') & fmask) ;
for ff=1:length(freqs)
  ind.(['gl' fstr{ff}])=find(lmask & fmask & p.band==freqs(ff));
end

% in BICEP2 A/B pairs are not adjacent in the channel order - scan the
% list looking for B partner to each A in each focal plane
[ind.a,ind.b] = pair_lookup(p, ind, find(amask));

% For backward compatibility, assign band==0 [dark] detectors to a particular
% frequency. Choose this by the median of the receiver.
band0mask = (p.band == 0 & amask);
band0a = cell(length(freqs), 1);
for ii=1:length(rxs)
  prxmask = (p.rx == rxs(ii));
  rxband0 = (band0mask & prxmask);
  if any(rxband0)
    % Get a representative frequency for this receiver (in reference to the
    % ordering in freqs).
    ff = find(median(p.band(prxmask)) == freqs);
    % Assign the band-0 detectors for this receiver
    band0a{ff} = [band0a{ff}; find(rxband0)];
  end
end
% Create frequency-specific A's and B's. Also merge in the band-0 detectors.
for ff=1:length(freqs)
  inda = intersect(ind.a, find(p.band == freqs(ff)));
  inda = unique([inda; band0a{ff}]);
  [inda,indb] = pair_lookup(p, ind, inda);
  ind.(['a' fstr{ff}]) = inda;
  ind.(['b' fstr{ff}]) = indb;
end

% find pairs where A and B are both good - we call these "rg" for
% "real good"
rg=find(p.flag(ind.a)==0&p.flag(ind.b)==0);
rg=[ind.a(rg);ind.b(rg)];
ind.rgl=intersect(ind.gl,rg);
ind.rgd=intersect(ind.gd,rg);
for ff=1:length(freqs)
  ind.(['rgl' fstr{ff}])=intersect(ind.(['gl' fstr{ff}]), rg);
end

% find A and B lists

% ind.gla and ind.glb are not a/b matched - they can be different lengths
ind.gla=intersect(ind.gl,ind.a);
ind.glb=intersect(ind.gl,ind.b);
for ff=1:length(freqs)
  ind.(['gl' fstr{ff} 'a'])=intersect(ind.(['gl' fstr{ff}]), ind.a);
  ind.(['gl' fstr{ff} 'b'])=intersect(ind.(['gl' fstr{ff}]), ind.b);
end

% these lists are a/b matched
[ind.la,dummy,inda]=intersect(ind.l,ind.a);
ind.lb=ind.b(inda);
[ind.rgla,dummy,inda]=intersect(ind.rgl,ind.a);
ind.rglb=ind.b(inda);
for ff=1:length(freqs)
  [lfreqa,dum,inda]=intersect(ind.(['l' fstr{ff}]), ind.a);
  ind.(['l' fstr{ff} 'a'])=lfreqa;
  ind.(['l' fstr{ff} 'b'])=ind.b(inda);
end
% Split this loop to maintain backward-compatibility with field order
for ff=1:length(freqs)
  [rglfreqa,dummy,inda]=intersect(ind.(['rgl' fstr{ff}]), ind.a);
  ind.(['rgl' fstr{ff} 'a']) = rglfreqa;
  ind.(['rgl' fstr{ff} 'b']) = ind.b(inda);
end

% Return row vectors!
fl=fieldnames(ind);
for i=1:length(fl)
  ind.(fl{i})=reshape(ind.(fl{i}),1,[]);
end

return

function [inda,indb] = pair_lookup(p,ind,inda)
  expt = get_experiment_name();

  switch expt
    % In BICEP3, the darks don't necessarily have a physical pair, so cull down
    % the list initially to only include the light detectors
    case 'bicep3'
      indafull = inda;
      inda = intersect(indafull, ind.l);
      indadark = setdiff(indafull, inda);
  end

  % in BICEP2 A/B pairs are not adjacent in the channel order - scan the
  % list looking for B partner to each A in each focal plane
  indb=zeros(size(inda));
  for i=1:length(inda)
    j=inda(i);
    pair=find(p.pix_phys_x==p.pix_phys_x(j) & p.pix_phys_y==p.pix_phys_y(j) ...
              & p.rx==p.rx(j), 2);
    indb(i)=pair(pair~=j);
  end

  switch expt
    % Now deal with the dark detectors. Since they don't have an obvious
    % physical pairing, we choose to pair the dark detectors in order. This
    % also maintains "pairs" on the same tile since there are 4 darks per
    % tile in BICEP3.
    case 'bicep3'
      inddark = zeros(length(ind.d)/2, 2);
      inddark(:,1) = intersect(ind.d, find(strcmp(p.pol, 'A')));
      inddark(:,2) = intersect(ind.d, find(strcmp(p.pol, 'B')));

      % Now match up the elements in indadark to inddark. Those rows are then
      % the pairs we need to append to each list.
      [dum,darkidx] = ismember(indadark, inddark(:,1));
      inda = [inda; inddark(darkidx,1)];
      indb = [indb; inddark(darkidx,2)];
      % Historically ind.a is always in GCP ascending order, so maintain
      % that, and fix-up ind.b to match.
      [inda,shufidx] = sort(inda);
      indb = indb(shufidx);
  end

  return
