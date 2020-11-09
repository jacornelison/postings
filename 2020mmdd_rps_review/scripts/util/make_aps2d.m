function [ad,aps2d,map]=make_aps2d(m,map,pure_b,pad,realimag,docross)
% [ad,aps2d,map]=make_aps2d(m,map,pure_b,pad,realimag,docross)
%
% pure_b - note if .QprojB etc are available will always get used
%
% pad=0 : don't pad maps
%    =1 : pad
%
% realimag=1 : plot left side imag(spec).*image(spec), right side real(spec).*real(spec)
%         =0 : plot real(spec.*conj(spec))
%
% docross - also calc the 2d cross spectra
%
% broken out from from reduc_makeaps2d

disp('calc 2d aps');

if ~exist('pure_b','var')||isempty(pure_b)
  pure_b='normal';
end
if ~exist('pad','var')||isempty(pad)
  pad=1;
end
if ~exist('realimag','var')||isempty(realimag)
  realimag=0;
end
if ~exist('docross','var')||isempty(docross)
  docross=1;
end

if(realimag == 1)
  display('Plotting imag on left and real on right')
end

% must follow reduc_makeaps to produce plots which represent what
% is going into the actual power spectra...

% pad the maps
if(pad)
  disp('padding maps before FT')
  for i=1:length(m)
    % this used to be one power higher than was used in reduc_makeaps
    % - now it is the same which means that new maps pager speclin plots
    % look a bit more pixelized versus old
    p=2^nextpow2(max([m(i).nx,m(i).ny]));
    [m(i),map(i,:)]=pad_map(m(i),map(i,:),p);
  end
end

% calc axis data
for i=1:length(m)
  ad(i)=calc_ad2([m(i).xdos,m(i).ydos],[m(i).nx,m(i).ny]);
end

for k=1:size(map,2)
  for j=1:size(map,1)
    
    % Q's included due to use in obscure maps pager plots
    % - naming of fields is a bit of a mess - should be TT, EE etc
    % but can't easily change

    % scale factor which is normally applied in calcspec
    sf=prod(ad(j).del_u);

    % do TT
    mw=map(j,k).Tw;
    ft=calc_map_fts(map(j,k),ad(j),mw,1,pure_b); % calc the weighted modes using map weight mw
    aps2d(j).T=specplot(ft.T,ft.T,sf,realimag);
    
    % do PP (EE, BB, EB)
    mw=map(j,k).Pw;
    ft=calc_map_fts(map(j,k),ad(j),mw,2,pure_b);
    aps2d(j).Q=specplot(ft.Q,ft.Q,sf,realimag);
    aps2d(j).U=specplot(ft.U,ft.U,sf,realimag);
    aps2d(j).E=specplot(ft.E,ft.E,sf,realimag);
    aps2d(j).B=specplot(ft.B,ft.B,sf,realimag);
    aps2d(j).EB=specplot(ft.E,ft.B,sf,realimag);
    % make matrix purified Q and U if they exist
    if isfield(ft,'QprojB')
      aps2d(j).QprojB=specplot(ft.QprojB,ft.QprojB,sf,realimag);
      aps2d(j).UprojB=specplot(ft.UprojB,ft.UprojB,sf,realimag);
    end

    % do TP (TE, TB)
    mw=gmean(map(j,k).Tw,map(j,k).Pw);
    ft=calc_map_fts(map(j,k),ad(j),mw,[1,2],pure_b);
    aps2d(j).TE=specplot(ft.T,ft.E,sf,realimag);
    aps2d(j).TB=specplot(ft.T,ft.B,sf,realimag);
    aps2d(j).TQ=specplot(ft.T,ft.Q,sf,realimag);
    aps2d(j).TU=specplot(ft.T,ft.U,sf,realimag);
  end

  if(docross)
    % append the inter row cross cross spectra in additional rows of aps array
    aps_offset = size(map,1);
    n = 1;
    for j=1:size(map,1)-1
      for c=j+1:size(map,1)
        
        % do only the base TT/EE/BB for now
        
        % do TT
        % make a common mask as gmean of two masks
        mw=gmean(map(j,k).Tw,map(c,k).Tw); 
        % mult each map by common mask and ft
        ft1=calc_map_fts(map(j,k),ad(j),mw,1,pure_b);
        ft2=calc_map_fts(map(c,k),ad(c),mw,1,pure_b);
        % if one ft smaller than other pad it (needed for high res LT)
        [ft1,ft2,ad(aps_offset+n)]=make_fts_same_size(ad(j),ad(c),ft1,ft2);
        % take the products of the two sets of modes
        aps2d(aps_offset+n).T=specplot(ft1.T,ft2.T,sf,realimag);
        
        % do PP
        mw=gmean(map(j,k).Pw,map(c,k).Pw);
        ft1=calc_map_fts(map(j,k),ad(j),mw,2,pure_b);
        ft2=calc_map_fts(map(c,k),ad(c),mw,2,pure_b);
        [ft1,ft2,ad(aps_offset+n)]=make_fts_same_size(ad(j),ad(c),ft1,ft2);
        aps2d(aps_offset+n).E=specplot(ft1.E,ft2.E,sf,realimag);
        aps2d(aps_offset+n).B=specplot(ft1.B,ft2.B,sf,realimag);
        %aps2d(aps_offset+n).EB=specplot(ft1.E,ft2.B,sf,realimag);
        %aps2d(aps_offset+n).BE=specplot(ft1.B,ft2.E,sf,realimag);
        
        % do TP
        %mw=gmean(map(j,k).Tw,map(c,k).Pw); 
        %ft1=calc_map_fts(map(j,k),ad,mw,1,pure_b);
        %ft2=calc_map_fts(map(c,k),ad,mw,2,pure_b);
        %aps2d(aps_offset+n).TE=specplot(ft1.T,ft2.E,sf,realimag);
        %aps2d(aps_offset+n).TB=specplot(ft1.T,ft2.B,sf,realimag);
        
        % do PT
        %mw=gmean(map(j,k).Pw,map(c,k).Tw);
        %ft1=calc_map_fts(map(j,k),ad,mw,2,pure_b);
        %ft2=calc_map_fts(map(c,k),ad,mw,1,pure_b);
        %aps2d(aps_offset+n).ET=specplot(ft1.E,ft2.T,sf,realimag);
        %aps2d(aps_offset+n).BT=specplot(ft1.B,ft2.T,sf,realimag);
        
        % increment to next slot
        n=n+1;
      end
    end
  end
end

if(0)
% mult to flat in l(l+1)C_l/2pi
% - are these now used anywhere?
ad.l_r=2*pi*ad.u_r;
sf=ad.l_r.*(ad.l_r+1)/(2*pi);
%disp('mult to flat in l(l+1)C_l/2pi');
for j=1:numel(map)
  aps2d(j).Tm=aps2d(j).T.*sf;
  if isfield(map(j),'Q')
    aps2d(j).Qm=aps2d(j).Q.*sf;
    aps2d(j).Um=aps2d(j).U.*sf;
    aps2d(j).Em=aps2d(j).E.*sf;
    aps2d(j).Bm=aps2d(j).B.*sf;
    aps2d(j).TQm=aps2d(j).TQ.*sf;
    aps2d(j).TUm=aps2d(j).TU.*sf;
    aps2d(j).TEm=aps2d(j).TE.*sf;
    aps2d(j).TBm=aps2d(j).TB.*sf;
    aps2d(j).EBm=aps2d(j).EB.*sf;
  end
  if isfield(map(j),'QprojB')
    aps2d(j).QprojBm=aps2d(j).Q.*sf;
    aps2d(j).UprojBm=aps2d(j).U.*sf;
  end
end
end

if(0)
if isfield(aps2d,'ET')
  n=1;
  for j=1:size(map,1)-1
    for c=j+1:size(map,1)
      aps2d(aps_offset+n).Tm=aps2d(aps_offset+n).T.*sf;
      aps2d(aps_offset+n).Qm=aps2d(aps_offset+n).Q.*sf;
      aps2d(aps_offset+n).Um=aps2d(aps_offset+n).U.*sf;
      aps2d(aps_offset+n).Em=aps2d(aps_offset+n).E.*sf;
      aps2d(aps_offset+n).Bm=aps2d(aps_offset+n).B.*sf;
      aps2d(aps_offset+n).TQm=aps2d(aps_offset+n).TQ.*sf;
      aps2d(aps_offset+n).TUm=aps2d(aps_offset+n).TU.*sf;
      aps2d(aps_offset+n).TEm=aps2d(aps_offset+n).TE.*sf;
      aps2d(aps_offset+n).TBm=aps2d(aps_offset+n).TB.*sf;
      aps2d(aps_offset+n).EBm=aps2d(aps_offset+n).EB.*sf;
      aps2d(aps_offset+n).QTm=aps2d(aps_offset+n).QT.*sf;
      aps2d(aps_offset+n).UTm=aps2d(aps_offset+n).UT.*sf;
      aps2d(aps_offset+n).ETm=aps2d(aps_offset+n).ET.*sf;
      aps2d(aps_offset+n).BTm=aps2d(aps_offset+n).BT.*sf;
      aps2d(aps_offset+n).BEm=aps2d(aps_offset+n).BE.*sf;
      n=n+1;
    end
  end
end
end

return

%%%%%%%%%%%%%%%%%%%%
function spec_out=specplot(spec_in_A,spec_in_B,sf,realimag)

if(~realimag)
  spec_out=sf*real(spec_in_A.*conj(spec_in_B));
else  
  width=size(spec_in_A,2);
  imagspec=sf*imag(spec_in_A).*imag(spec_in_B);
  realspec=sf*real(spec_in_A).*real(spec_in_B);
  spec_out(:,1:floor(width/2))=imagspec(:,1:floor(width/2));
  spec_out(:,floor(width/2)+1:width)=realspec(:,floor(width/2)+1:width);
end

return

