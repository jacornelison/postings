function fade_image(h,w)
% fade_image(h)
%
% fade an image away according to weight map m

cm=colormap; ncm=size(cm,1); ca=caxis;
d=get(h(1),'CData');
r=interp1(linspace(ca(1),ca(2),ncm),cm(:,1),d,[],'extrap');
g=interp1(linspace(ca(1),ca(2),ncm),cm(:,2),d,[],'extrap');
b=interp1(linspace(ca(1),ca(2),ncm),cm(:,3),d,[],'extrap');
im(:,:,1)=r; im(:,:,2)=g; im(:,:,3)=b;
w=repmat(w,[1,1,3]);
set(h(1),'CData',im.*w+(1-w));

return
