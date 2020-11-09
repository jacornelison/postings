function d=filter_scans2(d,fs,filttype,ind,p,difopt)
% d=filter_scans2(d,fs,filttype,ind, difopt)
%
% Filter half scans in the specified manner
%  this is taken from the original filter_scans and modified to allow for
%  regression of signal from diff pointing which causes leakage of T into P.
%
% Initially does not allow for any point source masking
%
% p is just the info from [p,ind]=get_array_info
% 
%  difopt is an option that controls which T map to get dif signal from
% (T, Tx, or Ty) and contains difopt.m, difopt.map and difopt.ukpervolt.
%  difopt also contains the fit coefficients for the regressor.
%  
% d=filter_scans(d,fs,'p3',ind.gla,p,difopt)
% d=filter_scans(d,fs,'p0',ind.glb,p,difopt)

  
if(~exist('ind'))
  ind=1:size(d.antenna0.bolo.mag,2);
end

if(~isfield(difopt,'type'))
  difopt.type='none';
end

if(~iscell(difopt.type))
  difopt.type={difopt.type};
end

switch filttype(1)
  case 'n'
    % do nothing
    
  case {'p', 'd', 's'}
    % subtract poly of required order
    d=polysub_scans(d,fs,str2num(filttype(2)),ind,p,difopt);

end

return

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d=polysub_scans(d,fs,porder,ind,p, difopt)
  
 %initialize dp
 nef=size(difopt.type,2);
 dp=zeros(size(d.ra,1),nef);
  
 %take data for that half scan for all channels considered
 v=double(d.antenna0.bolo.mag(:,ind));
       
  % loop over all channels
 for j=1:size(v,2)
  
   j
   % if even channel (pair-diff)
   if  mod(ind(j),2) == 0
   
     f = get_band_ind(p,ind,p.band(k));    
          
     % find ra/dec trajectory for this channel
     r=p.r(ind(j));
     th=p.theta(ind(j))-d.pointing.hor.dk;
     [y,x]=reckon(d.dec,d.ra,r,th+90);
     
     % Sample off the map at required coords
     for k=1:nef
     switch difopt.type{k}
      case 'T'
       %disp('T T T T T')
       sky=difopt.map(:,:,f).T;
       dp(:,k)=interp2nan(difopt.m.x_tic,difopt.m.y_tic,sky,x,y,'linear');
      case 'Tx'
       %disp('Tx Tx Tx Tx Tx')
       sky=difopt.map(:,:,f).Tx;
       dp(:,k)=interp2nan(difopt.m.x_tic,difopt.m.y_tic,sky,x,y,'linear');
      case 'Ty'
       %disp('Ty Ty Ty Ty Ty')
       sky=difopt.map(:,:,f).Ty;
       dp(:,k)=interp2nan(difopt.m.x_tic,difopt.m.y_tic,sky,x,y,'linear');
      case 'Te0'
       sky=difopt.map(:,:,f).Te0;
       dp(:,k)=interp2nan(difopt.m.x_tic,difopt.m.y_tic,sky,x,y,'linear');
     case 'Te1'
       sky=difopt.map(:,:,f).Te1;
       dp(:,k)=interp2nan(difopt.m.x_tic,difopt.m.y_tic,sky,x,y,'linear');
      case 'Te2'
       sky=difopt.map(:,:,f).Te2;
       dp(:,k)=interp2nan(difopt.m.x_tic,difopt.m.y_tic,sky,x,y,'linear');
     end
     end
     dp=dp/difopt.ukpervolt(f);
    
   end
   
   
   % for each half scan
   for i=1:length(fs.sf)
     
     s=fs.sf(i); e=fs.ef(i);scanind=s:e;
     azz=double(d.azoff(scanind));
     vh=v(scanind,j);
     
     % construct the regressor
     clear X
     for k=0:porder
       X(:,k+1)=azz.^k;
     end
     
     % if (type is T or Tx or Ty) and even channel (pair-diff)
     % add the diff regressors
     if ~strcmp(difopt.type{1},'none') &&  mod(ind(j),2) == 0 ...
           && mean(x(scanind)) < difopt.m.hx ... 
           && mean(x(scanind)) > difopt.m.lx
       
       X(:,porder+2:porder+1+nef)=dp(scanind,:);
       %plot(X(:,porder+2:porder+1+nef));pause
     end
     
     %compute baseline
     b=regress(vh,X);
     baseline=X*b;
     
     % subtract from unmasked data
     vh=vh-baseline;
     
     %plots
     if(0)
       clf;
       subplot(3,1,1)
       imagesc(difopt.m.x_tic, difopt.m.y_tic, sky) 
       subplot(3,1,1)
       hold on;
       plot(x,y,'.');
       title('WMAP T  and bolo trajectory')
       subplot(3,1,2)
       plot(vh);
       hold on
       plot(100*X(:,porder+2:end),'g');
       title('baseline sub bolo voltage(b),  100* WMAP template (g)')
       hold off
       subplot(3,1,3)
       plot(X(:,porder+2:end))
       title('diff pointing template')
       pause
     end
     
     end  % end loop over half scans
   
  % lscov can fit multiple columns in single shot
  % but annoyingly does not treat NaNs as missing values like regress
  % checked that one chan at a time regress above gives numerically
  % identical results
  %b=lscov(X,v);
  %v=v-X*b;

  % put poly sub data back in antenna0.bolo.mag
  d.antenna0.bolo.mag(scanind,ind(j))=vh;
  
end % end loop over channels

return

