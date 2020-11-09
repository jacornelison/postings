function out = rotate_3d(in, euler)
% out = rotate_3d(in, euler)
%
%   Rotate input 3-vectors according to Euler angles.
%
%   Euler angle convention is 1. rotate about z-axis by gamma
%                             2. rotate about y-axis by beta
%                             3. rotate about z-axis by alpha
%   
%   [Input arguments]
%     in      Input 3-vectors to be rotated. This argument can be
%             an array with shape [N,3], representing many vectors
%             to rotate.
%     euler   Euler angles for rotation transformation. This
%             argument can either have three elements, [alpha,
%             beta, gamma], in which case all input vectors will
%             undergo the same transformation. If this argument has
%             shape [N,3], then each input 3-vector will be rotated
%             by the corresponding set of Euler angles.
% 
%     euler(:,1) = alpha
%     euler(:,2) = beta
%     euler(:,3) = gamma
%
%   [Output]
%     out     Rotated 3-vectors. This array has the same shape as in.

% 2012-04-27 CAB

% Check size of argument in.
size_in = size(in);
if (ndims(in) > 2) | (size_in(end) ~= 3)
  disp('[rotate_3d] ERROR: argument in must have size [N,3]');
  return
end
nsamp = size_in(1);

% Check size of argument euler.
size_euler = size(euler);
if (ndims(euler) > 2) | (size_euler(end) ~= 3) | ...
      ((size_euler(1) ~= 1) & (size_euler(1) ~= nsamp))
  disp('[rotate_3d] ERROR: argument euler must have size [N,3] or [1,3]');
  return
end

% Pre-compute cos and sin of Euler angles.
dtor = pi / 180; % Convert degrees to radians.
ca = cos(dtor * euler(:,1));
sa = sin(dtor * euler(:,1));
cb = cos(dtor * euler(:,2));
sb = sin(dtor * euler(:,2));
cg = cos(dtor * euler(:,3));
sg = sin(dtor * euler(:,3));

% Create rotation matrices.
rot = zeros(nsamp, 3, 3);
rot(:,1,1) = ca .* cb .* cg - sa .* sg;
rot(:,1,2) = -1 * ca .* cb .* sg - sa .* cg;
rot(:,1,3) = ca .* sb;
rot(:,2,1) = sa .* cb .* cg + ca .* sg;
rot(:,2,2) = -1 * sa .* cb .* sg + ca .* cg;
rot(:,2,3) = sa .* sb;
rot(:,3,1) = -1 * sb .* cg;
rot(:,3,2) = sb .* sg;
rot(:,3,3) = cb;

% Apply rotation matrices.
if (nsamp > 1)
  out(:,1) = sum(squeeze(rot(:,1,:)) .* in, 2);
  out(:,2) = sum(squeeze(rot(:,2,:)) .* in, 2);
  out(:,3) = sum(squeeze(rot(:,3,:)) .* in, 2);
else
  % Hack to deal with weird behavior of the squeeze function.
  out(:,1) = sum(squeeze(rot(:,1,:))' .* in, 2);
  out(:,2) = sum(squeeze(rot(:,2,:))' .* in, 2);
  out(:,3) = sum(squeeze(rot(:,3,:))' .* in, 2);
end
