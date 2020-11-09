function d=make_utc_single_col(d)
  % d=make_utc_single_col(d);
  % utc comes back from reader as two column mjd/secinday - make
  % single col
  names=fieldnames(d);
  for i=1:length(names)
    if(isstruct(getfield(d,names{i})))
      d=setfield(d,names{i},make_utc_single_col(getfield(d,names{i})));
    else
      if(strfind(names{i},'utc')&size(getfield(d,names{i}),2)==2)
        tmp=getfield(d,names{i});
        mjddays=tmp(:,1);
        mjdtime=tmp(:,2);
        d=setfield(d,names{i},mjddays+mjdtime/86400);
      end
    end
  end

  return
