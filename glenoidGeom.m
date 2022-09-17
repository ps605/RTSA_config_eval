function scapula = glenoidGeom(R, hemi_gle_offsets, model_SSM, rhash)
%% Set up
% Load in and configure points of Scapula .stl
[x, y, z] = stlreadXYZ(['..\..\SSM\Scapulas\stl_aligned\' model_SSM '.stl']);

figure(10);

% Plot global coordinate system
x_hat=[0.1 0 0];
y_hat=[0 0.1 0];
z_hat=[0 0 0.1];

line([x_hat(1) 0],[x_hat(2) 0],[x_hat(3) 0], 'LineWidth',4,'Color','r'); % X - Red
line([y_hat(1) 0],[y_hat(2) 0],[y_hat(3) 0], 'LineWidth',4,'Color','y'); % Y - Yellow
line([z_hat(1) 0],[z_hat(2) 0],[z_hat(3) 0], 'LineWidth',4,'Color','g'); % Z - Green

% Plot .stl
patch(x,y,z,'b',...
    'FaceColor', [0.8 0.8 0.8],...
    'FaceAlpha', 0.75,...
    'EdgeColor', [0.65 0.65 0.65],...
    'EdgeAlpha', 0.5);

xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');

axis equal

view(3)
hold on;


%% Fit plane to glenoid rim points

% Glenoid rim points X-Y-Z
glenoid_points = importdata(['..\..\SSM\Scapulas\stl_aligned\' model_SSM 'rim_coords.txt'], ' ');
% Calculate the mean of the points
glenoid_barycentre = mean(glenoid_points);

% Get points to define coordinate system - Y-axis
[~, min_gl_p] = min(glenoid_points(:,2));
glenoid_inferior = glenoid_points(min_gl_p,:);

% Calculate vector between barrycentre and most inferior point
bc_to_inferior_p =  glenoid_inferior - glenoid_barycentre;

% Invert so it point superiorly and normalise
axial_normal = -(bc_to_inferior_p/norm(bc_to_inferior_p));

scatter3(glenoid_points(:,1), glenoid_points(:,2), glenoid_points(:,3), 'filled', 'o', 'cyan');
scatter3(glenoid_points(min_gl_p,1), glenoid_points(min_gl_p,2), glenoid_points(min_gl_p,3), 'filled', 'o', 'magenta');
scatter3(glenoid_barycentre(:,1), glenoid_barycentre(:,2), glenoid_barycentre(:,3), 'filled', 'o', 'magenta');

% Generate PointCloud of slected points and barycentre
glenoid_pointCloud = pointCloud(vertcat(glenoid_points, glenoid_barycentre));

figure (2);
pcshow(glenoid_pointCloud, 'MarkerSize', 100)

% Linear Regresion method to fit plane
x_gp = glenoid_points(:,1);
y_gp = glenoid_points(:,2);
z_gp = glenoid_points(:,3);

DM = [x_gp, y_gp, ones(size(z_gp))];
B = DM\z_gp;

% Create meshgrid of plane from Linear Regresion
[X,Y] = meshgrid(linspace(min(x_gp),max(x_gp),50), linspace(min(y_gp),max(y_gp),50));
Z = B(1)*X + B(2)*Y + B(3)*ones(size(X));

% Create point cloud Linear Regression plane (consistensy with following code)
plane_pointCloud = pointCloud([X(:), Y(:), Z(:)]);
% Fit plane to the Linear Regresion plane points
[glenoid_plane,~,~, ~] = pcfitplane(plane_pointCloud, 0.0001, 'MaxNumTrials', 1e6);

% Get normal to handle it later
% This is the Z-axis
glenoid_normal = glenoid_plane.Normal;

figure(10)

% Generate plane mesh and plot using Ax + By + Gz + D = 0
[plane_mesh_data.x_plane, plane_mesh_data.y_plane] = meshgrid(-0.1:0.01:0.1);
plane_mesh_data.z_plane = -1*(glenoid_plane.Parameters(1)*plane_mesh_data.x_plane ...
    + glenoid_plane.Parameters(2)*plane_mesh_data.y_plane ...
    + glenoid_plane.Parameters(4))/glenoid_plane.Parameters(3);

% Plot Plane
surf(plane_mesh_data.x_plane, plane_mesh_data.y_plane, plane_mesh_data.z_plane,...
    'FaceColor','g',...
    'FaceAlpha', 0.5,...
    'EdgeAlpha', 0.25)

%% Calculate scapular plane
scap_pointCloud = pointCloud([x(:), y(:), z(:)]);

% Linear Regression method to fit plane
x_sp = x(:);
y_sp = y(:);
z_sp = z(:);

DM = [x_sp, y_sp, ones(size(z_sp))];
B = DM\z_sp;

% Create meshgrid of plane from Linear Regresion
[X,Y] = meshgrid(linspace(min(x_sp),max(x_sp),50), linspace(min(y_sp),max(y_sp),50));
Z = B(1)*X + B(2)*Y + B(3)*ones(size(X));

% Create point cloud Linear Regression plane (consistensy with following code)
scap_plane_pointCloud = pointCloud([X(:), Y(:), Z(:)]);
% Fit plane to the Linear Regresion plane points
[scap_plane,~,~, ~] = pcfitplane(scap_plane_pointCloud, 0.0001, 'MaxNumTrials', 1e6);

% Generate plane mesh and plot using Ax + By + Gz + D = 0
[plane_mesh_data.x_plane, plane_mesh_data.y_plane] = meshgrid(-0.1:0.01:0.1);
plane_mesh_data.z_plane = -1*(scap_plane.Parameters(1)*plane_mesh_data.x_plane ...
    + scap_plane.Parameters(2)*plane_mesh_data.y_plane ...
    + scap_plane.Parameters(4))/scap_plane.Parameters(3);

figure;
pcshow(scap_pointCloud, 'MarkerSize',20);
hold on;
surf(plane_mesh_data.x_plane, plane_mesh_data.y_plane, plane_mesh_data.z_plane,...
    'FaceColor','y',...
    'FaceAlpha', 0.5,...
    'EdgeAlpha', 0.25)
%% Project Points onto glenoid plane (minimisation problem)
% "Most inferior" glenoid point onto glenoid plane
plane_parameters = glenoid_plane.Parameters;
% Constraint function (plane)
f_con = @(y_m)plane_func_d(y_m, plane_parameters);
% Cost function (distance)
J_inferior = @(y_m)sqrt((glenoid_inferior(1) - y_m(1))^2 + (glenoid_inferior(2) - y_m(2))^2 + (glenoid_inferior(3) - y_m(3))^2);
% Initial Condition (point on plane)
y_m_0 = [plane_mesh_data.x_plane(1) plane_mesh_data.y_plane(1) plane_mesh_data.z_plane(1)];

options = optimset('MaxIter', 100, 'TolFun', 1e-4);
% Run fmincon
[glenoid_inferior_on_plane, ~] = fmincon(J_inferior,...
    y_m_0,...
    [],...
    [],...
    [],...
    [],...
    [],...
    [],...
    f_con,...
    options);

scatter3(glenoid_inferior_on_plane(1), glenoid_inferior_on_plane(2), glenoid_inferior_on_plane(3), 'filled','o','magenta');

% Glenoid barycentre onto glenoid plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This updates the value of glenoid_barycentre to the projection (x, y, z)
% Cost function (distance)
J_barycentre = @(y_m)sqrt((glenoid_barycentre(1) - y_m(1))^2 + (glenoid_barycentre(2) - y_m(2))^2 + (glenoid_barycentre(3) - y_m(3))^2);
% Run fmincon
[glenoid_barycentre, ~] = fmincon(J_barycentre,...
    y_m_0,...
    [],...
    [],...
    [],...
    [],...
    [],...
    [],...
    f_con,...
    options);

scatter3(glenoid_barycentre(1), glenoid_barycentre(2), glenoid_barycentre(3), 'filled','o','magenta');

% Check if Barycentre sits on plane. Should be -> 0
bary_plane = glenoid_plane.Parameters(1)*glenoid_barycentre(1) +...
    glenoid_plane.Parameters(2)*glenoid_barycentre(2) +...
    glenoid_plane.Parameters(3)*glenoid_barycentre(3) + ...
    glenoid_plane.Parameters(4);

if bary_plane >= 1e-4
    disp('Error: You fucked up')
    keyboard
end


%% Check which way the Palne norm (Z axis) is pointing

% Get the barycentre of glenoid verteces
x_sca_mean = mean(x,'all');
y_sca_mean = mean(y,'all');
z_sca_mean = mean(z,'all');

scap_barycentre = [x_sca_mean y_sca_mean z_sca_mean];

scatter3(scap_barycentre(1), scap_barycentre(2), scap_barycentre(3),'cyan', 'filled','o')

% Pass .stl x-y-z to structure to then export
stl_scap = struct('x', x,...
    'y', y,...
    'z', z);

% What is the angle between the norm and the vector from the norm origine to glenoid vertex barycentre

% Vector from glenoid barycentre to humeral verteces barrycentre
gle_bary_vec = scap_barycentre - glenoid_barycentre;

% Connect with line to visualise vector
line([glenoid_barycentre(1) scap_barycentre(1)],...
    [glenoid_barycentre(2) scap_barycentre(2)],...
    [glenoid_barycentre(3) scap_barycentre(3)], ...
    'LineWidth',4,'Color','cyan');

norm_check_angle = vrrotvec(gle_bary_vec, glenoid_normal);

% If angle between vectors is 0<gonia<=90 flip normal about plane
% This is the Z-axis
if norm_check_angle(4) >= 0 && norm_check_angle(4) < pi/2

    % Negate Z normal to flip about plane to point out of the glenoid
    glenoid_normal = - glenoid_normal;
    % Project point normal from plane to distance R along Z normal
    glenoid_plane_normals.z_p = (glenoid_barycentre + R.*glenoid_normal);
    scatter3(glenoid_plane_normals.z_p(1), glenoid_plane_normals.z_p(2), glenoid_plane_normals.z_p(3),'green', 'filled','o','MarkerEdgeColor','black')

    % Connect with line to visualise normal and projection
    line([glenoid_barycentre(1) glenoid_plane_normals.z_p(1)],...
        [glenoid_barycentre(2) glenoid_plane_normals.z_p(2)],...
        [glenoid_barycentre(3) glenoid_plane_normals.z_p(3)], ...
        'LineWidth',4,'Color','green');

else
    % Project point normal from plane to distance R along Z normal
    glenoid_plane_normals.z_p = (glenoid_barycentre + R.*glenoid_normal);
    scatter3(glenoid_plane_normals.z_p(1), glenoid_plane_normals.z_p(2), glenoid_plane_normals.z_p(3),'green', 'filled','o','MarkerEdgeColor','black')

    % Connect with line to visualise normal and projection
    line([glenoid_barycentre(1) glenoid_plane_normals.z_p(1)],...
        [glenoid_barycentre(2) glenoid_plane_normals.z_p(2)],...
        [glenoid_barycentre(3) glenoid_plane_normals.z_p(3)], ...
        'LineWidth',4,'Color','green');

end

% Plot Barycetntre where cup will be placed
scatter3(glenoid_barycentre(1), glenoid_barycentre(2), glenoid_barycentre(3), 'black','filled','o')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if plane norm and glenoid_inferior_on_plane are perpendicular
bc_to_inf_glen = (glenoid_inferior_on_plane - glenoid_barycentre)/norm(glenoid_inferior_on_plane - glenoid_barycentre);

if dot(bc_to_inf_glen, glenoid_normal) < 1e-10
    % Negate so that axis normal is pointing superiorly
    glenoid_plane_y_n = -bc_to_inf_glen;
elseif dot(bc_to_inf_glen, glenoid_normal) >= 1e-10
    disp(' Error: Glenoid plane Y and Z axes not perpendicular (dot(bc_to_inf_glen, glenoid_normal) >= 1e-10)');
    keyboard
end
%% Handle hemisphere at glenoid plane

% Define depth/geometry of the "glenoid head" (phi<theta/4 for anything smaller than
% hemisphere)
theta = (0:0.01:1)*2*pi;
phi = (0:0.01:1)*pi/2;

% Need to rototransalte it so the the vertex is sitting on the point of
% interest
[THETA,PHI]=meshgrid(theta,phi);
X1=R.*cos(THETA).*sin(PHI) + glenoid_barycentre(1);
Y1=R.*sin(THETA).*sin(PHI) + glenoid_barycentre(2);
Z1=R.*cos(PHI) + glenoid_barycentre(3);

figure(10);
hemisphere_gle = surf(X1,Y1,Z1,...
    'FaceColor',[ 1 1 0],...
    'FaceAlpha', 0.75,...
    'EdgeColor', [0 0 0 ],...
    'EdgeAlpha', 0.1);
axis equal


% Quck workaround to rotate - Use the graphic object handler and then
% extract the point data X-Y-Z

% Find axis and angle of rotation between plane normal and where hemisphere
% is ploted about Z-axis [0 0 1]
glen_rot = vrrotvec([0 0 1], glenoid_normal);
% Rotate hemisphere about plane normal axis to aligned on plane out of
% plane normal

rotate(hemisphere_gle,glen_rot(1:3),rad2deg(glen_rot(4)), glenoid_barycentre)

% This is the Glenoid plane Y-axis
glenoid_plane_normals.y_p = glenoid_barycentre + glenoid_plane_y_n(1:3)*R;
scatter3(glenoid_plane_normals.y_p(1),glenoid_plane_normals.y_p(2), glenoid_plane_normals.y_p(3),'yellow','filled','o','MarkerEdgeColor','black')

line([glenoid_barycentre(1) glenoid_plane_normals.y_p(1)],...
    [glenoid_barycentre(2) glenoid_plane_normals.y_p(2)],...
    [glenoid_barycentre(3) glenoid_plane_normals.y_p(3)], ...
    'LineWidth',4,'Color','yellow');

%% Extract baseline parameters and data points for cup, plane and humerus
% Extract relevant data from hemisphere cup placed on resection barycentre
% and from there make positional changes as below.

glenoid_plane_normals.y_n = glenoid_plane_y_n; % Superior/inferior
glenoid_plane_normals.z_n = glenoid_normal; % Out of plane
glenoid_plane_normals.x_n = -cross(glenoid_plane_normals.z_n,glenoid_plane_normals.y_n); % Anterior/Posterior

% Plot axis of final norm from barycentre
glenoid_plane_normals.x_p = glenoid_barycentre + glenoid_plane_normals.x_n(1:3)*R;
scatter3(glenoid_plane_normals.x_p(1),glenoid_plane_normals.x_p(2), glenoid_plane_normals.x_p(3),'red','filled','o','MarkerEdgeColor','black')

line([glenoid_barycentre(1) glenoid_plane_normals.x_p(1)],...
    [glenoid_barycentre(2) glenoid_plane_normals.x_p(2)],...
    [glenoid_barycentre(3) glenoid_plane_normals.x_p(3)], ...
    'LineWidth',4,'Color','red');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check they're mutually perpendicular
if dot(glenoid_plane_normals.x_n, glenoid_plane_normals.y_n ) > 1e-10 && ...
        dot(glenoid_plane_normals.y_n, glenoid_plane_normals.z_n) > 1e-10

    disp('Error: Plane normals not perpendicular')

    keyboard

end

% All displacements are defined on the glenoid plane now based on the
% variable: glenoid_plane_normals


%% Change position of the cup

% 1) Position on resection surface (superior/inferior, anterior/posterior)
% 2) Offset from resection surface
% 3) Radius (wrt to glenosphere for congruent fit)
% 4) Depth
% 5) Version/Inclination

% Get mesh data for the hemisphere from Visualisation Object. This will
% then be continuesly updated through "hemisphere" Surface variable

hemi_gle_mesh_data.X = hemisphere_gle.XData;
hemi_gle_mesh_data.Y = hemisphere_gle.YData;
hemi_gle_mesh_data.Z = hemisphere_gle.ZData;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The vector that will be rotated needs to be a column vector so it is
% transposed here (') then back again after rotations are finished
% Translate cup centre to originate from origin
CoR_glen = glenoid_barycentre;

%% Version/Inclination of cup about hemisphere normal axes

% Antero-/Postero- version (about Proximal/Distal axis)
rotate(hemisphere_gle,...
    glenoid_plane_normals.y_n,...
    hemi_gle_offsets.y_ant_retro_version,...
    glenoid_barycentre)

% Rotate CoR_glen about glenoid_plane_normals.y_n
R_y = axang2rotm([glenoid_plane_normals.y_n deg2rad(hemi_gle_offsets.y_ant_retro_version)]);

% Need to rotate glenoid_plane_normals.x_n axis after first rotation
glenoid_plane_normals.x_n_r1 = R_y*glenoid_plane_normals.x_n';
glenoid_plane_normals.x_n_r1 = glenoid_plane_normals.x_n_r1';

% Supero-/Infero- inclination (about Anterior/Posterior axis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For this rotation it will be rotated about the new X axis after first Z
% rotation to keep topological meaning for the cup orientation. Can be
% thought of as the orientation of the cup as it sat on the resection plane
% and then rotated about first axes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Negate (-) hemi_gle_offsets.x_sup_inf_incl so that positive is superior
hemi_gle_offsets.x_sup_inf_incl = - hemi_gle_offsets.x_sup_inf_incl;

rotate(hemisphere_gle,...
    glenoid_plane_normals.x_n_r1,...
    hemi_gle_offsets.x_sup_inf_incl,...
    glenoid_barycentre)

R_x = axang2rotm([glenoid_plane_normals.x_n_r1 deg2rad(hemi_gle_offsets.x_sup_inf_incl)]);

% Rotate glenoid_plane_normals.y_n axis after second rotation
glenoid_plane_normals.y_n_r1 = R_x*glenoid_plane_normals.y_n';
glenoid_plane_normals.y_n_r1 = glenoid_plane_normals.y_n_r1';


% Get transformed axes offsets from origin 
glenoid_plane_normals.theta(1) = atan2(norm(cross(glenoid_plane_normals.x_n_r1,[1 0 0])),dot(glenoid_plane_normals.x_n_r1,[1 0 0]));
glenoid_plane_normals.theta(2) = atan2(norm(cross(glenoid_plane_normals.y_n_r1,[0 1 0])),dot(glenoid_plane_normals.y_n_r1,[0 1 0]));
glenoid_plane_normals.theta(3) = 0;
%% Position on glenoid surface (anterior/posterior, base offset, superior/inferior)

% X - Anterior / Posterior offsets
hemisphere_gle.XData = hemisphere_gle.XData + glenoid_plane_normals.x_n(1)*hemi_gle_offsets.x_ant_post;
hemisphere_gle.YData = hemisphere_gle.YData + glenoid_plane_normals.x_n(2)*hemi_gle_offsets.x_ant_post;
hemisphere_gle.ZData = hemisphere_gle.ZData + glenoid_plane_normals.x_n(3)*hemi_gle_offsets.x_ant_post;

% Adjust barrycentre to now be Joint CoR
CoR_glen = CoR_glen + glenoid_plane_normals.x_n*hemi_gle_offsets.x_ant_post;

% Y - Proximal / Distal offsets
hemisphere_gle.XData = hemisphere_gle.XData + glenoid_plane_normals.y_n(1)*hemi_gle_offsets.y_prox_dist;
hemisphere_gle.YData = hemisphere_gle.YData + glenoid_plane_normals.y_n(2)*hemi_gle_offsets.y_prox_dist;
hemisphere_gle.ZData = hemisphere_gle.ZData + glenoid_plane_normals.y_n(3)*hemi_gle_offsets.y_prox_dist;

% Adjust barrycentre to now be Joint CoR
CoR_glen = CoR_glen + glenoid_plane_normals.y_n*hemi_gle_offsets.y_prox_dist;

% Z - Base offset
hemisphere_gle.XData = hemisphere_gle.XData + glenoid_plane_normals.z_n(1)*hemi_gle_offsets.z_base_off;
hemisphere_gle.YData = hemisphere_gle.YData + glenoid_plane_normals.z_n(2)*hemi_gle_offsets.z_base_off;
hemisphere_gle.ZData = hemisphere_gle.ZData + glenoid_plane_normals.z_n(3)*hemi_gle_offsets.z_base_off;

% Adjust barrycentre to now be Joint CoR
CoR_glen = CoR_glen + glenoid_plane_normals.z_n*hemi_gle_offsets.z_base_off;

scatter3(CoR_glen(1),CoR_glen(2), CoR_glen(3),'magenta','filled','o','MarkerEdgeColor','black')


% % % keyboard

%% Create scapula/glenoid structure to output for manipulation

scapula = struct('glenoid_plane_normals', glenoid_plane_normals,... % x-y-z normal vectors and end-points
    'hemi_gle_mesh_data', hemi_gle_mesh_data,...                    % DEFAULT data for glenoid hemisphere on glenoid plane before any rototranslation
    'hemi_gle_offsets', hemi_gle_offsets,...                        % glenoid hemisphere rototranslation offsets
    'glenoid_barycentre', glenoid_barycentre,...                    % centre of DEFAULT hemisphere (CoR of joint effectivly - needs the hemi_gle_offsets translation only to be applied along normals)
    'hemisphere_gle', hemisphere_gle,...                            % rototranslated glenoid hemisphere
    'glenoid_plane', glenoid_plane,...                              % glenoid plane parameters
    'plane_mesh_data', plane_mesh_data,...
    'R', R,...                                                      % hemisphere radius (m)
    'CoR_glen', CoR_glen,...                                        % adjusted barycentre now as CoR
    'stl_scap', stl_scap);

stlwrite_user(['..\..\OpenSim\In\Geometry\gle_' rhash '.stl'],...
    hemisphere_gle.XData,...
    hemisphere_gle.YData,...
    hemisphere_gle.ZData,...
    'mode','ascii',...
    'triangulation','f');


end
