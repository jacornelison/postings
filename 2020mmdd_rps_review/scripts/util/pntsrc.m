function V=pntsrc(par,u,v)
% V=pntsrc(par,u,v)
%
% Calculate the complex visibilities expected for a source
% location in par and a set of baselines u,v
%
% par(1) = flux
% par(2:3) = l,m position on field

% Dot product of src location and baseline vectors
a=2*pi*(par(2)*u+par(3)*v);

% Expected visibility
V=par(1)*complex(cos(a),sin(a));
