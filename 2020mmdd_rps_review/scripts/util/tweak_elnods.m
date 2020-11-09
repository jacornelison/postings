function enc=tweak_elnods(d,en,sampratio)
% enc=tweak_elnods(d,en,sampratio)
%
% Tweak el nod start/end indices to center on the el nod. The el nod
% region identified by the feature bits contains extra padding at the
% start and is slightly truncated at the end. This function identi-
% fies the el nod based on actual mount motion and shifts the start/
% end indices to center on it, without changing the length.

enc=en;
el=d.antenna0.tracker.actual(:,2);
loop_closed=bitand(d.antenna0.pmac.mtr_stat,2^5)==0;
nsnap=d.array.frame.nsnap;

% Mask of times NaNed out at edges of deconvolution --
% Don't want to run into these
dcmask=false(sampratio*length(loop_closed));
% Loop over MCEs as defined by get_xferfunc
for i=1:length(d.tf(:))
  dm=d.tf(i).dcmask;
  % Check if dm is empty due to unused MCE
  if ~isempty(dm)
    for j=1:length(dm.sf)
      dcmask(dm.sf(j):dm.ef(j))=true;
    end
  end
end
% Turn into slow mask
dcmask=any(reshape(dcmask,sampratio,[]),1);

for i=1:length(en.s)
  % Grab the elevation value for the marked el nod period
  eltmp=el(en.s(i):en.e(i));
  idx=(1:length(eltmp)) - 1;

  % Find 'combine' number for this el nod
  nsnap_good=nanmedian(nsnap(en.s(i):en.e(i)));

  % If the el nod has too few samples, it's garbage and following
  % code will fail. Bail out here and keep initial s/e indices. 
  if length(idx)<4
    continue
  end

  % If the servo loop wasn't closed, something strange was going on.
  if any(~loop_closed(en.s(i):en.e(i)))
    continue
  end

  % Identify the samples that are in the +/- peaks by percentile.
  % This is more robust than picking out a single min/max sample.
  % But if this fails, fall back on min/max.
  pcnt=percentile(eltmp,[0.25,0.75]);
  if diff(pcnt)==0 || sum(eltmp<pcnt(1))==0 || sum(eltmp>pcnt(2))==0
    pcnt=percentile(eltmp,[0 1]);
  end
  % El nod center is avg sample # of + and - peak regions.
  idx_cent=(mean(idx(eltmp<=pcnt(1))) + mean(idx(eltmp>=pcnt(2)))) / 2;

  % Calculate back to el nod start and finish.
  n=en.e(i)-en.s(i);
  enc.s(i)=round(idx_cent+en.s(i)-n/2);
  enc.e(i)=enc.s(i)+n;

  % If we've gone off the reservation, keep the old indices.
  if (enc.s(i)<1) || (enc.e(i)>length(el)) || (enc.e(i)<=enc.s(i))
    enc.s(i)=en.s(i);
    enc.e(i)=en.e(i);
    continue
  end

  % Make sure we don't have any changes in sample rate.
  % If first/last frame has different nsnap, shrink
  % adjusted el nod region until it doesn't.
  % Don't shrink past original boundaries (in that case,
  % it's not our fault).
  while(enc.s(i)<en.s(i) && (nsnap(enc.s(i))~=nsnap_good || dcmask(enc.s(i))))
    enc.s(i)=enc.s(i)+1;
  end
  while(enc.e(i)>en.e(i) && (nsnap(enc.e(i))~=nsnap_good || dcmask(enc.e(i))))
    enc.e(i)=enc.e(i)-1;
  end
end

% Introduce sf,ef - start,end fast sample
% (In BICEP 1st fast sample is contemporaneous with 1st slow sample
% according to antenna0.time.utcfast/slow)
enc.sf=(enc.s-1)*sampratio+1;
enc.ef=(enc.e-1)*sampratio+1;

return

