function R=make_rotations(alpha,beta,gamma)  
  % First euler rotation: rotate normal vector by azimuthal angle alpha:
% Ra=[cosd(alpha)  sind(alpha) 0; ...
%     -sind(alpha) cosd(alpha) 0; ...
%     0            0           1];
Ra=zeros(3,3,length(alpha));
Ra(1,1,:)=cosd(alpha);
Ra(1,2,:)=sind(alpha);
Ra(2,1,:)=-sind(alpha);
Ra(2,2,:)=cosd(alpha);
Ra(3,3,:)=1;
  
% Second euler rotation: rotate normal vector by declination beta
% Rb=[cosd(beta) 0 sind(beta); ...
%     0          1 0 ; ...
%    -sind(beta) 0 cosd(beta)];
Rb=zeros(3,3,length(alpha));
Rb(1,1,:)=cosd(beta);
Rb(1,3,:)=sind(beta);
Rb(2,2,:)=1;
Rb(3,1,:)=-sind(beta);
Rb(3,3,:)=cosd(beta);
  
% Third euler rotation: rotate normal vector by "roll" of the mirror angle
% gamma:
% Rg=[cosd(gamma)  sind(gamma) 0; ...
%     -sind(gamma) cosd(gamma) 0; ...
%     0           0          1];
Rg=zeros(3,3,length(alpha));
Rg(1,1,:)=cosd(gamma);
Rg(1,2,:)=sind(gamma);
Rg(2,1,:)=-sind(gamma);
Rg(2,2,:)=cosd(gamma);
Rg(3,3,:)=1;  

R=multiprod(Rg,multiprod(Rb,Ra));
