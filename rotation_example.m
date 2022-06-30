%Create the meshgrid
Xinit = 0:3125:800000;
Yinit = 0:3125:800000;
[Xq,Yq] = meshgrid(Xinit,Yinit);
Z = sin(Xq/40000) .* cos(Yq/40000);                     % ‘Z’ To Provide A Surface
figure(11)
mesh(Xq,Yq,Z) 
hold on
% Original
grid on
XY = [Xq(:) Yq(:)];                                     % Create Matrix Of Vectors
theta=-180; %TO ROTATE CLOCKWISE BY X DEGREES
R=[cosd(theta) -sind(theta); sind(theta) cosd(theta)]; %CREATE THE MATRIX
rotXY=XY*R'; %MULTIPLY VECTORS BY THE ROT MATRIX 
Xqr = reshape(rotXY(:,1), size(Xq,1), []);
Yqr = reshape(rotXY(:,2), size(Yq,1), []);
%SHIFTING
Xqrs = Xqr+175.5362;
Yqrs = Yqr+175.5362;
figure(11)
mesh(Xqrs, Yqrs, Z)                                     % Rotated & Shifted
grid on