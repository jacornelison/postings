function [T,R]=qu2linpol(Q,U,ad)
% [T,R]=qu2linpol(Q,U)
%
% Make pol vectors from Q,U maps
% If no return args given instead plots as quiver style plot
%
% T = polarization angle
% R = polarization magnitude

T=0.5*atan2(U,Q);
R=sqrt(U.^2+Q.^2);

if(nargout==0)
  [u,v]=pol2cart(T,R);
  
  % Make sub-sample array
  ind=1:10:ad.N_pix;
  x=ad.t_val_deg{1}(ind);

  quiver(x,x,+u(ind,ind),+v(ind,ind),0.5,'g.');
  hold on;
  quiver(x,x,-u(ind,ind),-v(ind,ind),0.5,'g.');
  hold off;
  axis square; axis tight;
end
