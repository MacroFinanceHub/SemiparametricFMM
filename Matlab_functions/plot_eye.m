%%%%%%%%%%%%%%%%%%
% function to plot spherical data as a heatmap on a circle
% effectively does an orthographic azimuthal projection with center at 
% latitude=0 (north pole)
%
% Since it is orthographic, the perspective is from an infinite distance
% away, which causes the region near the equator to be more compressed
% relative to the poles, which are more spread out.
%
% Note that the most one should plot on the same plot is a hemisphere,
% since this leads to a unique, invertible projection.  Another plot could
% be done for the southern hempisphere of complete spherical data
%
%
% zz=intensities on spherical surface to plot, vector of length t1*t2
% (defined below)
% radius = constant for radius of sphere
% latitude = set of latitudes (actually zenith angles since on (0,pi) not (-pi/2, pi/2)
%           corresponding to data zz, length t1
% longitude = set of longitudes (azimuthal angles) corresponding to data zz, length t2
% zlim=vector of [minz, maxz] defining range of data for heatmap -- i.e.
% zz>maxz will be pure red, zz<minz will be pure blue, and then the numbers
% in between follow a grid defined by the color wheel.
%
% min_radius=minimum radius to plot -- white out everything inside.
% 
% Sometimes the mesh works, sometimes it doesn't -- seemingly randomly it
% will give garbage, then work fine, then give garbage again!!!! 4/3/13
%  To solve this problem, do:
% set(gcf, 'Renderer', 'zbuffer');
%
% Also, note that heatmap plots with negative y values up instead of down,
% so reverse y values from polar coordinates to make picture (but not axes)
% right

function plot_eye(zz,radius,latitude,longitude,zlim,min_radius)

if (nargin<5)
    zlim=[min(zz),max(zz)];
end;
if (nargin<6)
    min_radius=2000;  %%% default for eye data 04/13
end;

T_long=length(longitude);
T_lat=length(latitude);
latitude=repmat(latitude,1,T_long);
longitude=repmat(longitude',T_lat,1);
x=radius*sin(latitude).*cos(longitude);
y=-radius*sin(latitude).*sin(longitude);

xx=reshape(x,1,T_long*T_lat);
yy=reshape(y,1,T_long*T_lat);

xlin=linspace(min(xx),max(xx),1000);
ylin=linspace(min(yy),max(yy),1000);
[XX,YY]=meshgrid(xlin,ylin);

ZZ3=griddata(xx,yy,zz,XX,YY,'cubic');
ZZ2=(XX.^2+YY.^2)>min_radius^2;
ZZ=nan(size(XX));
ZZ(ZZ2==1)=ZZ3(ZZ2==1);
mesh(XX,YY,ZZ)
caxis(zlim)
view(0,270)
%imagesc(ZZ)
  %% Why the heck does mesh not work now???????????????????????????
