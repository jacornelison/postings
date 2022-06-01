function [cbin, xbin, ybin, phasebin, Mpoints, local_phase, Tmap] = moon_illumination(mpos,spos,ang_diam,np_ang,ppd,datadir)

% Calculate Illumination per pixel
% We'll have Ra/Dec/HA/angular diameter

% Create a sphere with some pixels, each with normals pointing outward
% from the center.

%%

O = [0;0;0];
X = [1;0;0];
Y = [0;1;0];
Z = [0;0;1];

% point cloud
samples = 20000;
Mpoints = -fibonacci_sphere(samples,0);
samples = size(Mpoints,2);
R = rodmatrix(np_ang,Z);
Mpoints_latlon = R*Mpoints;
Mlat = asind(Mpoints_latlon(2,:));
Mlat(abs(Mlat)>64) = 64;
mvec = [0;0;-1];
svec = Z;
R = rodmatrix(np_ang,-Z);
np_vec = R*Y;
eq_vec = R*X;

% Rotate into RA/DEC
% R1 = rodmatrix(mpos(2)*0,X);
% R2 = rodmatrix(mpos(1)*0,Y);
% Mpoints = R2*R1*Mpoints;

xpos = spos(2)-mpos(2);
ypos = spos(1)-mpos(1);
R1 = rodmatrix(xpos,[-1;0;0]);
R2 = rodmatrix(ypos,Y);
svec = R2*R1*svec;

STO = acosd(dot(mvec,svec));



if 0
    
    % Z X Y
fig = figure(2);
clf; hold on;
plot3(-Mpoints(3,:)*0.5,Mpoints(1,:)*0.5,Mpoints(2,:)*0.5,'.')
quiver3(0,0,0,1,0,0,'Color','k')
quiver3(0,0,0,0,1,0,'Color','r')
quiver3(0,0,0,0,0,1,'Color','b')

%quiver3(0,0,0,svec(1),svec(2),svec(3))
%quiver3(0,0,0,np_vec(3),np_vec(1),np_vec(2),'Color','m')
%quiver3(0,0,0,eq_vec(3),eq_vec(1),eq_vec(2),'Color','m')
ang = 0;
V = Z*0.5;
V2 = V;
for i = -90:90%1:360
    V = [V rodmatrix(i,np_vec)*Z*0.5];
    V2 = [V2 rodmatrix(i,eq_vec)*Z*0.5];
end

plot3(V(3,:),V(1,:),V(2,:),'k.','MarkerSize',3)    
plot3(V2(3,:),V2(1,:),V2(2,:),'k.','MarkerSize',3)

grid on;
xlabel('Z')
ylabel('X')
zlabel('Y')
view(105,10)
axis equal
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])
end
% sun-moon direction
%dirvec = svec-mvec*0.0025;
dirvec = svec;
dirvec = dirvec/norm(dirvec);

%%
% angle in the moon-sun plane
% Compute the local solar phase based on the Sun's position WRT to Moon.
% local_phase = atan2(dot(Mpoints,repmat(svec,1,size(Mpoints,2))),...
%     dot(Mpoints,repmat(mvec,1,size(Mpoints,2))))*180/pi;
local_phase = dot(Mpoints,repmat(dirvec,1,size(Mpoints,2)));

% Flat projection toward the observer
saz = atan2(dot(svec,Z),dot(svec,eq_vec))*180/pi;
Maz = atan2(dot(Mpoints,repmat(Z,1,samples)),dot(Mpoints,repmat(eq_vec,1,samples)))*180/pi;
local_phase = saz-Maz+180;


% Create an Temperature profile based on Bruce et al. 1965 Fig 34/35.
load(fullfile(datadir,'moon_temp_data.mat'));
Tmap = griddata(phase,lat,BT,wrapTo180(local_phase),Mlat,'v4');


if 0
       
fig = figure(3);
clf; hold on;
scatter(phase,lat,5,BT,'filled')
colorbar()
colormap default
grid on
scatter(local_phase,Mlat,15,Tmap,'filled')

fig = figure(4);
clf; hold on;
ind = isnan(Tmap);
scatter(Mpoints(1,ind),Mpoints(2,ind),14,local_phase(ind),'filled')
%scatter(Mpoints(1,:),Mpoints(2,:),14,Mlat,'filled')
set(gca,'xdir','reverse')
colorbar()

end



xbin = (-0.3:1/ppd/2:0.3);
ybin = (-0.3:1/ppd/2:0.3);
phasebin = grid_map_simple(Mpoints(1,:)*ang_diam/2,Mpoints(2,:)*ang_diam/2,local_phase,xbin,ybin);

%cbin = grid_map_simple(Mpoints(1,:)*ang_diam/2,Mpoints(2,:)*ang_diam/2,local_phase,xbin,ybin);
cbin = grid_map_simple(Mpoints(1,:)*ang_diam/2,Mpoints(2,:)*ang_diam/2,Tmap,xbin,ybin);
%cbin = grid_map_simple(Mpoints(1,:)*ang_diam/2,Mpoints(2,:)*ang_diam/2,asind(Mpoints_latlon(1,:)),xbin,ybin);


function M = fibonacci_sphere(samples,clip)

samples = 2*samples;
inds = 1:samples;

phi = acos(1-2*inds/samples);
theta = pi*(1+sqrt(5))*inds;

x = cos(theta).*sin(phi);
y = sin(theta).*sin(phi);
z = cos(phi);

ind = z>clip/90;
M = [x(ind); y(ind); z(ind);];


function zbin = grid_map_simple(x,y,z,xbin,ybin)

inrange = @(A,B,C) B<=A & A<=C;
zbin = NaN(length(xbin),length(ybin));

dx = median(diff(xbin));
dy = median(diff(ybin));
for xi = 1:length(xbin)
    for yi = 1:length(ybin)
        ind = find(inrange(x,xbin(xi)-dx,xbin(xi)+dx) &...
            inrange(y,ybin(yi)-dy,ybin(yi)+dy));
        if ~isempty(ind)
            
                zbin(xi,yi) = nanmean(z(ind));

        end
    end
end










