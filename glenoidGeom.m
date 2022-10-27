function scapula = glenoidGeom(R, hemi_gle_offsets, model_SSM, rhash, flag_correctVersion, flag_correctInclination, flag_correct12mm)
%% Set up
% Load in and configure points of Scapula .stl
% NOTE: not the most efficient way of handling the .stl
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
glenoid_points = importdata(['..\..\SSM\Scapulas\stl_aligned\' model_SSM '_rim_coords.txt'], ' ');
% Calculate the mean of the points
glenoid_barycentre = mean(glenoid_points);

% Get points to define coordinate system - Y-axis
D = pdist(glenoid_points);
D = squareform(D);
[~,I] = max(D(:));
[max_point_idx_1, max_point_idx_2] = ind2sub(size(D),I);
max_point_idx = [max_point_idx_1, max_point_idx_2];

% % % [~, min_gl_p] = min(glenoid_points(:,2));
% % % glenoid_inferior = glenoid_points(min_gl_p,:);
% % % 
% % % % Calculate vector between barrycentre and most inferior point
% % % bc_to_inferior_p =  glenoid_inferior - glenoid_barycentre;
% % % 
% % % % Invert so it point superiorly and normalise
% % % axial_normal = -(bc_to_inferior_p/norm(bc_to_inferior_p));

scatter3(glenoid_points(:,1), glenoid_points(:,2), glenoid_points(:,3), 'filled', 'o', 'cyan');
% scatter3(glenoid_points(min_gl_p,1), glenoid_points(min_gl_p,2), glenoid_points(min_gl_p,3), 'filled', 'o', 'magenta');
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
[gle_plane_mesh_data.x_plane, gle_plane_mesh_data.y_plane] = meshgrid(-0.1:0.01:0.1);
gle_plane_mesh_data.z_plane = -1*(glenoid_plane.Parameters(1)*gle_plane_mesh_data.x_plane ...
    + glenoid_plane.Parameters(2)*gle_plane_mesh_data.y_plane ...
    + glenoid_plane.Parameters(4))/glenoid_plane.Parameters(3);

% Plot Plane
surf(gle_plane_mesh_data.x_plane, gle_plane_mesh_data.y_plane, gle_plane_mesh_data.z_plane,...
    'FaceColor','g',...
    'FaceAlpha', 0.25,...
    'EdgeAlpha', 0)

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
[sca_plane_mesh_data.x_plane, sca_plane_mesh_data.y_plane] = meshgrid(-0.1:0.01:0.1);
sca_plane_mesh_data.z_plane = -1*(scap_plane.Parameters(1)*sca_plane_mesh_data.x_plane ...
    + scap_plane.Parameters(2)*sca_plane_mesh_data.y_plane ...
    + scap_plane.Parameters(4))/scap_plane.Parameters(3);

% figure;
% pcshow(scap_pointCloud, 'MarkerSize',20);
% hold on;
surf(sca_plane_mesh_data.x_plane, sca_plane_mesh_data.y_plane, sca_plane_mesh_data.z_plane,...
    'FaceColor','y',...
    'FaceAlpha', 0.25,...
    'EdgeAlpha', 0)

%% Calculate supraspinatus fossa base vector

load('fossa_base.mat');
principal_cmp = pca([x(fossa_base.vertices(:)), y(fossa_base.vertices(:)), z(fossa_base.vertices(:))]);
fossa_vector = principal_cmp(:,1)';
% Plot 1st principle component vector
%fossa_point_i = [x(fossa_base.vertices(9)), y(fossa_base.vertices(9)), z(fossa_base.vertices(9))];
% fossa_point_f = fossa_point_i + fossa_vector.*0.1;
fossa_point_f = glenoid_barycentre + fossa_vector.*R;



% line([fossa_point_i(1) fossa_point_f(1)], [fossa_point_i(2) fossa_point_f(2)], [fossa_point_i(3) fossa_point_f(3)], 'LineWidth',4,'Color','cyan'); 
line([glenoid_barycentre(1) fossa_point_f(1)], [glenoid_barycentre(2) fossa_point_f(2)], [glenoid_barycentre(3) fossa_point_f(3)], 'LineWidth',4,'Color','cyan'); 
scatter3(fossa_point_f(1), fossa_point_f(2), fossa_point_f(3), 'filled', 'cyan', 'o','MarkerEdgeColor','black')


%% Project Points onto glenoid plane (minimisation problem)
% "Most inferior" glenoid point onto glenoid plane
plane_parameters = glenoid_plane.Parameters;
% Constraint function (plane)
f_con = @(y_m)plane_func_d(y_m, plane_parameters);

for i_max_point = 1:2
    % Cost function (distance)
    J_inferior = @(y_m)sqrt((glenoid_points(max_point_idx(i_max_point),1) - y_m(1))^2 + (glenoid_points(max_point_idx(i_max_point),2) - y_m(2))^2 + (glenoid_points(max_point_idx(i_max_point),3) - y_m(3))^2);
    % Initial Condition (point on plane)
    y_m_0 = [gle_plane_mesh_data.x_plane(1) gle_plane_mesh_data.y_plane(1) gle_plane_mesh_data.z_plane(1)];

    options = optimset('MaxIter', 100, 'TolFun', 1e-4);
    % Run fmincon
    [max_point_plane, ~] = fmincon(J_inferior,...
        y_m_0,...
        [],...
        [],...
        [],...
        [],...
        [],...
        [],...
        f_con,...
        options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is currently note used as now using plane interscection but keep for
% ease if wanted to re-introduce
%    glenoid_max_points(i_max_point,:) = max_point_plane;
    scatter3(max_point_plane(1), max_point_plane(2), max_point_plane(3), 'filled','o','magenta');
end
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

%% Calculate vector of glenoid and scapula plane intersection

[~, intersect_v] = plane_intersect(glenoid_plane.Parameters(1:3),...
                             [gle_plane_mesh_data.x_plane(1) gle_plane_mesh_data.y_plane(1) gle_plane_mesh_data.z_plane(1)],...
                             scap_plane.Parameters(1:3),...
                             [sca_plane_mesh_data.x_plane(1) sca_plane_mesh_data.y_plane(1) sca_plane_mesh_data.z_plane(1)]);

% Normalise intersection vector
intersect_v = intersect_v/norm(intersect_v);

% Project point from barycentre along intersect axis for visualisation
pl_p = glenoid_barycentre + R*intersect_v;

% % % Connect with line to visualise normal and projection
% %     line([glenoid_barycentre(1) pl_p(1)],...
% %         [glenoid_barycentre(2) pl_p(2)],...
% %         [glenoid_barycentre(3) pl_p(3)], ...
% %         'LineWidth',4,'Color','green');

% Check for orientation of vector with respect to Y axis
intersect_v_angle = vrrotvec([0 1 0], intersect_v);
intersect_v_angle_deg = rad2deg(intersect_v_angle(4));

%% Check which way the Plane norm (Z axis) is pointing

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
%vec_max_glen_points = (glenoid_max_points(1,:) - glenoid_max_points(2,:))/norm(glenoid_max_points(1,:) - glenoid_max_points(2,:));

if dot(intersect_v, glenoid_normal) < 1e-10
    % NEED to make sure axis is pointingg superiorly
    if intersect_v_angle_deg < 90 && intersect_v_angle_deg > - 90
        glenoid_plane_y_n = intersect_v;
    else
        glenoid_plane_y_n = -intersect_v;
    end
elseif dot(vec_max_glen_points, glenoid_normal) >= 1e-10
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
glenoid_plane_normals.y_p = glenoid_barycentre + R.*glenoid_plane_y_n(1:3);
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

%% Calculate version and inclination correction angles from fossa vector angle; And 12 mm rule from most inferior point on glenoid Y-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% YZ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fossa_vector onto glenoid YZ plane (glenoid_plane_normals.x_n)
% Constraint function (plane)
delta = -(glenoid_plane_normals.x_n(1)*glenoid_plane_normals.z_p(1) + glenoid_plane_normals.x_n(2)*glenoid_plane_normals.z_p(2) + glenoid_plane_normals.x_n(3)*glenoid_plane_normals.z_p(3));
plane_parameters = [glenoid_plane_normals.x_n delta];
f_con = @(y_m)plane_func_d(y_m, plane_parameters);
y_m_0 = glenoid_barycentre;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This updates the value of glenoid_barycentre to the projection (x, y, z)
% Cost function (distance)
J_fossa_point_f = @(y_m)sqrt((fossa_point_f(1) - y_m(1))^2 + (fossa_point_f(2) - y_m(2))^2 + (fossa_point_f(3) - y_m(3))^2);
% Run fmincon
[fossa_point_f_YZ, ~] = fmincon(J_fossa_point_f,...
    y_m_0,...
    [],...
    [],...
    [],...
    [],...
    [],...
    [],...
    f_con,...
    options);

% Calculate Inclination correction angle
 
% Calculate unit vector from barycentre to projected point
fossa_correction_v.YZ = (fossa_point_f_YZ - glenoid_barycentre)/norm(fossa_point_f_YZ - glenoid_barycentre);
% Calcuate angle between correction angle about x_n
fossa_correction_ang.YZ = vrrotvec(glenoid_plane_normals.z_n, fossa_correction_v.YZ);
% Push correction vector point to R from barycentre
fossa_point_f_YZ = glenoid_barycentre + fossa_correction_v.YZ*R;

% Visualise Inclination angle
scatter3(fossa_point_f_YZ(1), fossa_point_f_YZ(2), fossa_point_f_YZ(3), 'filled','o','cyan', 'MarkerEdgeColor','black');
version_poly = [glenoid_barycentre; fossa_point_f_YZ; glenoid_plane_normals.z_p];
patch(version_poly(:,1), version_poly(:,2) , version_poly(:,3), 'r');
line([glenoid_barycentre(1) fossa_point_f_YZ(1)],...
    [glenoid_barycentre(2) fossa_point_f_YZ(2)],...
    [glenoid_barycentre(3) fossa_point_f_YZ(3)], ...
    'LineWidth',4,'Color','g');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% XZ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fossa_vector onto glenoid XZ plane (glenoid_plane_normals.y_n)
% Constraint function (plane)
delta = -(glenoid_plane_normals.y_n(1)*glenoid_plane_normals.z_p(1) + glenoid_plane_normals.y_n(2)*glenoid_plane_normals.z_p(2) + glenoid_plane_normals.y_n(3)*glenoid_plane_normals.z_p(3));
plane_parameters = [glenoid_plane_normals.y_n delta];
f_con = @(y_m)plane_func_d(y_m, plane_parameters);
y_m_0 = glenoid_barycentre;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This updates the value of glenoid_barycentre to the projection (x, y, z)
% Cost function (distance)
J_fossa_point_f = @(y_m)sqrt((fossa_point_f(1) - y_m(1))^2 + (fossa_point_f(2) - y_m(2))^2 + (fossa_point_f(3) - y_m(3))^2);
% Run fmincon
[fossa_point_f_XZ, ~] = fmincon(J_fossa_point_f,...
    y_m_0,...
    [],...
    [],...
    [],...
    [],...
    [],...
    [],...
    f_con,...
    options);

% Calculate Version correction angle
 
% Calculate unit vector from barycentre to projected point
fossa_correction_v.XZ = (fossa_point_f_XZ - glenoid_barycentre)/norm(fossa_point_f_XZ - glenoid_barycentre);
% Calcuate angle between correction angle about x_n
fossa_correction_ang.XZ = vrrotvec(glenoid_plane_normals.z_n, fossa_correction_v.XZ);
% Push correction vector point to R from barycentre
fossa_point_f_XZ = glenoid_barycentre + fossa_correction_v.XZ*R;

% Visualise Version angle
scatter3(fossa_point_f_XZ(1), fossa_point_f_XZ(2), fossa_point_f_XZ(3), 'filled','o','cyan', 'MarkerEdgeColor','black');
version_poly = [glenoid_barycentre; fossa_point_f_XZ; glenoid_plane_normals.z_p];
patch(version_poly(:,1), version_poly(:,2) , version_poly(:,3), 'y');
line([glenoid_barycentre(1) fossa_point_f_XZ(1)],...
    [glenoid_barycentre(2) fossa_point_f_XZ(2)],...
    [glenoid_barycentre(3) fossa_point_f_XZ(3)], ...
    'LineWidth',4,'Color','g');

%%%%%%%%%%%%%%%%%%%%%%% Clean up correction angles %%%%%%%%%%%%%%%%%%%%%%%%
correction_angles.x_sup_inf_incl        = rad2deg(fossa_correction_ang.YZ(4));
correction_angles.y_ant_retro_version   = rad2deg(fossa_correction_ang.XZ(4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 12 mm rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project a point inferiorly along glenoid -ive Y-axis
p_point = glenoid_barycentre-glenoid_plane_normals.y_n*R*2;
% Get smallest Euclidian distance of glenoid rim points from projected
% point on -ive Y-axis. Not exact bur close ennough
min_rim_points = vecnorm((glenoid_points - p_point), 2 , 2);
[~, inf_point_idx] = min(min_rim_points);

% Calculate distances
inf_point = glenoid_points(inf_point_idx, :);
d_inferior = norm(inf_point - glenoid_barycentre);

correction_displacement.y_prox_dist = d_inferior - 0.012;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1st Rotation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Supero-/Infero- inclination (about Anterior/Posterior axis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For this rotation it will be rotated about the new X axis after first Z
% rotation to keep topological meaning for the cup orientation. Can be
% thought of as the orientation of the cup as it sat on the resection plane
% and then rotated about first axes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_correctInclination == true
    % Set sup_inf_inclination to correction value calculated
    hemi_gle_offsets.x_sup_inf_incl = correction_angles.x_sup_inf_incl;
else
    % Negate (-) hemi_gle_offsets.x_sup_inf_incl so that positive is superior
    hemi_gle_offsets.x_sup_inf_incl = - hemi_gle_offsets.x_sup_inf_incl;
end

rotate(hemisphere_gle,...
    glenoid_plane_normals.x_n,...
    hemi_gle_offsets.x_sup_inf_incl,...
    glenoid_barycentre)

% Rotation matrix
R_x = axang2rotm([glenoid_plane_normals.x_n deg2rad(hemi_gle_offsets.x_sup_inf_incl)]);

% Rotate glenoid_plane_normals.y_n axis after second rotation
glenoid_plane_normals.y_n_r1 = R_x*glenoid_plane_normals.y_n';
glenoid_plane_normals.y_n_r1 = glenoid_plane_normals.y_n_r1';

glenoid_plane_normals.z_n_r1 = R_x*glenoid_plane_normals.z_n';
glenoid_plane_normals.z_n_r1 = glenoid_plane_normals.z_n_r1';

ppy = glenoid_barycentre + R*glenoid_plane_normals.y_n_r1;
scatter3(ppy(1), ppy(2), ppy(3), 'yellow', 'filled');

ppz = glenoid_barycentre + R*glenoid_plane_normals.z_n_r1;
scatter3(ppz(1), ppz(2), ppz(3), 'green', 'filled');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2nd Rotation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_correctVersion == true
    % Set ant_retro_version to correction value calculated
    hemi_gle_offsets.y_ant_retro_version = correction_angles.y_ant_retro_version;
end

% Antero-/Postero- version (about Proximal/Distal axis)
rotate(hemisphere_gle,...
    glenoid_plane_normals.y_n_r1,...
    hemi_gle_offsets.y_ant_retro_version,...
    glenoid_barycentre)

% Rotation matrix
R_y = axang2rotm([glenoid_plane_normals.y_n_r1 deg2rad(hemi_gle_offsets.y_ant_retro_version)]);
R_z = axang2rotm([glenoid_plane_normals.z_n_r1 deg2rad(hemi_gle_offsets.y_ant_retro_version)]);

% Need to rotate glenoid_plane_normals.x_n axis after first rotation
glenoid_plane_normals.x_n_r1 = R_y*glenoid_plane_normals.x_n';
glenoid_plane_normals.x_n_r1 = glenoid_plane_normals.x_n_r1';

glenoid_plane_normals.z_n_r2 = R_y*glenoid_plane_normals.z_n_r1';
glenoid_plane_normals.z_n_r2 = glenoid_plane_normals.z_n_r2';

ppx = glenoid_barycentre + R*glenoid_plane_normals.x_n_r1;
scatter3(ppx(1), ppx(2), ppx(3), 'red', 'filled');

ppz = glenoid_barycentre + R*glenoid_plane_normals.z_n_r2;
scatter3(ppz(1), ppz(2), ppz(3), 'green', 'filled');


% Get transformed axes orientation offsets from origin 

% Final rotation matrix of glenosphere axes
RM = [glenoid_plane_normals.x_n_r1;...
    glenoid_plane_normals.y_n_r1;...
    glenoid_plane_normals.z_n_r2];

% Calculate eular angles of fianl RM in global
ZYX_Euler_ang = rotm2eul(RM, 'ZYX');

% Invert calcutated angles
glenoid_plane_normals.theta(1) = - ZYX_Euler_ang(3);
glenoid_plane_normals.theta(2) = - ZYX_Euler_ang(2);
glenoid_plane_normals.theta(3) = - ZYX_Euler_ang(1);
%% Position on glenoid surface (anterior/posterior, base offset, superior/inferior)

% X - Anterior / Posterior offsets
hemisphere_gle.XData = hemisphere_gle.XData + glenoid_plane_normals.x_n(1)*hemi_gle_offsets.x_ant_post;
hemisphere_gle.YData = hemisphere_gle.YData + glenoid_plane_normals.x_n(2)*hemi_gle_offsets.x_ant_post;
hemisphere_gle.ZData = hemisphere_gle.ZData + glenoid_plane_normals.x_n(3)*hemi_gle_offsets.x_ant_post;

% Adjust barrycentre to now be Joint CoR
CoR_glen = CoR_glen + glenoid_plane_normals.x_n*hemi_gle_offsets.x_ant_post;

% Y - Proximal / Distal offsets
if flag_correct12mm == true
    hemi_gle_offsets.y_prox_dist = - correction_displacement.y_prox_dist;
end

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

ppy = CoR_glen + R*glenoid_plane_normals.y_n_r1;
scatter3(ppy(1), ppy(2), ppy(3), 'yellow', 'filled');

ppz = CoR_glen + R*glenoid_plane_normals.z_n_r2;
scatter3(ppz(1), ppz(2), ppz(3), 'green', 'filled');
% % % keyboard

%% Create scapula/glenoid structure to output for manipulation

scapula = struct('glenoid_plane_normals', glenoid_plane_normals,... % x-y-z normal vectors and end-points
    'hemi_gle_mesh_data', hemi_gle_mesh_data,...                    % DEFAULT data for glenoid hemisphere on glenoid plane before any rototranslation
    'hemi_gle_offsets', hemi_gle_offsets,...                        % glenoid hemisphere rototranslation offsets
    'glenoid_barycentre', glenoid_barycentre,...                    % centre of DEFAULT hemisphere (CoR of joint effectivly - needs the hemi_gle_offsets translation only to be applied along normals)
    'hemisphere_gle', hemisphere_gle,...                            % rototranslated glenoid hemisphere
    'glenoid_plane', glenoid_plane,...                              % glenoid plane parameters
    'plane_mesh_data', sca_plane_mesh_data,...
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
