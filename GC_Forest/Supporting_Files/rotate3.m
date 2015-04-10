%% Rotates 3D surface based on input angles

function [X,Y,Z] = rotate3(M,xdeg,ydeg,zdeg)
M=M*[1 0 0; 0 cosd(xdeg) sind(xdeg); 0 -sind(xdeg) cosd(xdeg)];
M=M*[cosd(ydeg) 0 -sind(ydeg); 0 1 0; sind(ydeg) 0 cosd(ydeg)];
M=M*[cosd(zdeg) sind(zdeg) 0; -sind(zdeg) cosd(zdeg) 0; 0 0 1];
X=M(:,1);
Y=M(:,2);
Z=M(:,3);
end
