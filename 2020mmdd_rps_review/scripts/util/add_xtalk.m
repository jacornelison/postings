function d=add_xtalk(d,fs,en,xl,p,ind,do_relgain)
% d=add_xtalk(d,fs,en,xl,p,ind,do_relgain)
%
% Adds inductive crosstalk to the TODs. 
%
% d,fs,en  = regular data structures
% xl       = level of crosstalk; 0 = (no xtalk), i.e. .01, add 1 percent of a channel's
%            signal to its mux_row neighbors 
% p,ind    = regular get_array_info structures
% do_xtalk = 1 to include ADU/airmass correction (default)
%            0 to do xtalk level in units of source map signal
%
% CDS 12/7/12

if ~xl
  % return without doing anything
  return
end
if ~exist('do_relgain','var') || isempty(do_relgain)
  do_relgain=true;
end

% adjust the size of xl.  can be full per-det or 1 number for the array
switch length(xl)
 case 1
  % one value per array.  expand to full size
  xl=xl.*ones(size(ind.e));
 case length(ind.e)
  % one value per det. don't do anything
 otherwise
  error('dont know how to expand xl (crosstalk levels)')
end

% Inductive crosstalk 

% According to RWA:
% Inductive xtalk just goes by mux row number, so channel A will see its channel B
% polarization pair, and the upstream channel B.  Channel B will see its channel A
% polarization buddy and the downstream channel A. There is a wraparound at the end
% of the mux row.

% do this per-rx for keck
for rx=unique(p.rx)'

  % strip down the ind for future loop
  indrx=strip_ind(ind,find(p.rx==rx));

  mux_col=unique(p.mce_col);
  for k=1:numel(mux_col)
  
    % pull out this column, in which pairs cross talk among each other
    colind=ind.e(p.mce_col==mux_col(k) & p.rx==rx);
    v=d.mce0.data.fb(:,colind);
  
    % store a copy of these timestreams
    v2=v;

    if do_relgain
      % find the gains (ADU/airmass) for these channels
      g=nanmean(en.g(:,colind,2),1);
  
      % set NaN gains to median gain, just to avoid headaches
      g(~isfinite(g))=nanmedian(g);
    else
      g=ones(size(colind));
    end
  
    % these are the mux row numbers corresponding to this block (single column) of
    % timestreams 
    mux_row=p.mce_row(colind);
    rowsize=numel(mux_row);

    for j=1:rowsize
      % This is the timestream we will add to its neighbors
      mux_row_0=mux_row(j);
      % This is its upstream neighbor
      mux_row_up=mod(mux_row_0+1,rowsize);
      % This is its downstream neighbor
      mux_row_down=mod(mux_row_0-1,rowsize);
    
      % indices within the column that has been pulled out
      ind0=find(mux_row==mux_row_0);
      indup=find(mux_row==mux_row_up);
      inddown=find(mux_row==mux_row_down);

      % index within the ind.e structure
      indxl=find(p.mce_row==mux_row_0 & p.mce_col==mux_col(k) & p.rx==rx);
    
      % add in xtalk
      % timestreams will keep the same mean in the map region, though it won't matter for
      % mean subtracted maps or for signal only sims.
      % Multiply by the ADU/airmass gains before adding in crosstalk, because this is
      % what happens inside the MCE, then undo by dividing by the gain
      mapind=make_mapind(d,fs);
      v0=v(:,ind0);
      vu=v(:,indup)-nanmean(v(mapind,indup));
      vd=v(:,inddown)-nanmean(v(mapind,inddown));

      v2(:,ind0)=(vu*g(indup)*xl(indxl) + vd*g(inddown)*xl(indxl) + v0*g(ind0))/g(ind0);
    
    end

    % Copy this crosstalk added mux column back into the data structure
    d.mce0.data.fb(:,colind)=v2;

  end % column
end % rx

return



    
