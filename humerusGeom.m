function humerus = humerusGeom(R, hemi_cup_offsets,rhash)
%% Set up
% Load in and configure points of Humerurs .stl
[x, y, z] = stlreadXYZ('..\..\OpenSim\In\Geometry\humerus_resected_manifold_closed.stl');

figure(1);

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
    'FaceAlpha', 0.5,...
    'EdgeColor', [0.65 0.65 0.65],...
    'EdgeAlpha', 0.5);

xlabel('X-axis');
ylabel('Y-axis');
zlabel('z-axis');

axis equal

view(3)
hold on;

% % % % % Pause to select points on .stl where "resectin" was made and load in from
% % % % % plot
% % % % keyboard;
% % % %
% % % % resection_points =  evalin('base', 'cursor_info');
p = load('resection_points.mat');

%% Fit plane to resection points

% Calculate the mean of the points
resection_barycentre = mean(vertcat(p.resection_points(:).Position));

% Generate PointCloud of slected points and barycentre
resection_pointCloud = pointCloud(vertcat(p.resection_points(:).Position, resection_barycentre));

% figure (2);
% pcshow(resection_pointCloud, 'MarkerSize', 100)

% Fit plane to the points
[resection_plane,~,~, ~] = pcfitplane(resection_pointCloud, 0.001);
% Get normal to handle it later
resection_normal = resection_plane.Normal;

figure(1)

% Generate plane mesh and plot using Ax + By + Gz + D = 0
[plane_mesh_data.x_plane, plane_mesh_data.y_plane] = meshgrid(-0.1:0.01:0.1);
plane_mesh_data.z_plane = -1*(resection_plane.Parameters(1)*plane_mesh_data.x_plane ...
    + resection_plane.Parameters(2)*plane_mesh_data.y_plane ...
    + resection_plane.Parameters(4))/resection_plane.Parameters(3);

% Plot Plane
surf(plane_mesh_data.x_plane, plane_mesh_data.y_plane, plane_mesh_data.z_plane,...
    'FaceColor','g',...
    'FaceAlpha', 0.5,...
    'EdgeAlpha', 0.25)


%% Handle hemisphere at humeral resection

% Define depth/geometry of the "cup" (phi<theta/4 for anything smaller than
% hemisphere)
theta = (0:0.01:1)*2*pi;
phi = (0:0.01:1)*pi/4;


% Plot Barycetntre where cup will be placed
scatter3(resection_barycentre(1), resection_barycentre(2), resection_barycentre(3), 'cyan','filled','o')

% Project point normal to the plane from Barycentre
% Check what direction it points
resection_plane_normals.y_p = (resection_barycentre + R.*resection_normal);
scatter3(resection_plane_normals.y_p(1), resection_plane_normals.y_p(2), resection_plane_normals.y_p(3),'yellow', 'filled','o','MarkerEdgeColor','black')
% Connect with line to visualise normal and projection
line([resection_barycentre(1) resection_plane_normals.y_p(1)],...
    [resection_barycentre(2) resection_plane_normals.y_p(2)],...
    [resection_barycentre(3) resection_plane_normals.y_p(3)], ...
    'LineWidth',4,'Color','yellow');


%% Check which way the norm is pointing

% Get the barycentre of humerus verteces
x_hum_mean = mean(x,'all');
y_hum_mean = mean(y,'all');
z_hum_mean = mean(z,'all');

hum_bary = [x_hum_mean y_hum_mean z_hum_mean];

scatter3(hum_bary(1), hum_bary(2), hum_bary(3),'black', 'filled','o')

% Pass .stl x-y-z to structure to then export
stl_hum = struct('x', x,...
    'y', y,...
    'z', z);


% What is the angle between the norm and the vector from the norm origine to humeral vertex barycentre

% Vector from resection barycentre to humeral verteces barrycentre
hum_bary_vec = hum_bary - resection_barycentre;

% Connect with line to visualise vector
line([resection_barycentre(1) hum_bary(1)],...
    [resection_barycentre(2) hum_bary(2)],...
    [resection_barycentre(3) hum_bary(3)], ...
    'LineWidth',4,'Color','cyan');

norm_check_angle = vrrotvec(hum_bary_vec, resection_normal);

% If angle between vectors is 0<gonia<=90 flip normal abou plane
if norm_check_angle(4) >= 0 && norm_check_angle(4) < pi/2
    resection_normal = - resection_normal;

    line([resection_barycentre(1) resection_plane_normals.y_p(1)],...
        [resection_barycentre(2) resection_plane_normals.y_p(2)],...
        [resection_barycentre(3) resection_plane_normals.y_p(3)], ...
        'LineWidth',4,'Color','yellow');
end

%% continue
% Check if Barycentre sits on plane. Should be -> 0
bary_plane = resection_plane.Parameters(1)*resection_barycentre(1) +...
    resection_plane.Parameters(2)*resection_barycentre(2) +...
    resection_plane.Parameters(3)*resection_barycentre(3) + ...
    resection_plane.Parameters(4);

if bary_plane >= 1e-4
    disp('Barycentre not on plane')
    keyboard
end

% Need to rototransalte it so the the vertex is sitting on the point of
% interest
[THETA,PHI]=meshgrid(theta,phi);
X1=R.*cos(THETA).*sin(PHI) + resection_barycentre(1);
Y1=R.*sin(THETA).*sin(PHI) + resection_barycentre(2);
Z1=R.*cos(PHI) + resection_barycentre(3);

figure(1);
hemisphere_hum = surf(X1,Y1,Z1,...
    'FaceColor',[ 1 1 0],...
    'FaceAlpha', 0.75,...
    'EdgeColor', [0 0 0 ],...
    'EdgeAlpha', 0.1);
axis equal



% Quck workaround to rotate - Use the graphic object handler and then
% extract the point data X-Y-Z

% Find axis and angle of rotation between plane normal and where hemisphere
% is ploted about Z-axis [0 0 1]

hemi_rot = vrrotvec([0 0 1], -resection_normal);
resection_plane_normals.x_p = resection_barycentre + hemi_rot(1:3)*R;
scatter3(resection_plane_normals.x_p(1), resection_plane_normals.x_p(2), resection_plane_normals.x_p(3),'red','filled','o','MarkerEdgeColor','black')

line([resection_barycentre(1) resection_plane_normals.x_p(1)],...
    [resection_barycentre(2) resection_plane_normals.x_p(2)],...
    [resection_barycentre(3) resection_plane_normals.x_p(3)], ...
    'LineWidth',4,'Color','red');

% Rotate hemisphere about plane normal axis to aligned on plane out of
% plane normal

rotate(hemisphere_hum,hemi_rot(1:3),rad2deg(hemi_rot(4)), resection_barycentre)

rotX = hemisphere_hum.XData;
rotY = hemisphere_hum.YData;
rotZ = hemisphere_hum.ZData;

hemisphere_hum.XData = rotX+resection_normal(1)*R;
hemisphere_hum.YData = rotY+resection_normal(2)*R;
hemisphere_hum.ZData = rotZ+resection_normal(3)*R;


clear rotX rotY rotZ
%% Extract baseline parameters and data points for cup, plane and humerus
% Extract relevant data from hemisphere cup placed on resection barycentre
% and from there make positional changes as below.

resection_plane_normals.x_n = hemi_rot(1:3); % Anterior/Posterior
resection_plane_normals.y_n = resection_normal; % Out of plane
resection_plane_normals.z_n = cross(resection_plane_normals.x_n,resection_plane_normals.y_n); %Proximal/Distal

% Plot axis of final norm from barycentre
resection_plane_normals.z_p = resection_barycentre + resection_plane_normals.z_n(1:3)*R;
scatter3(resection_plane_normals.z_p(1),resection_plane_normals.z_p(2), resection_plane_normals.z_p(3),'green','filled','o','MarkerEdgeColor','black')

line([resection_barycentre(1) resection_plane_normals.z_p(1)],...
    [resection_barycentre(2) resection_plane_normals.z_p(2)],...
    [resection_barycentre(3) resection_plane_normals.z_p(3)], ...
    'LineWidth',4,'Color','green');

% check they're mutually perpendicular
if dot(resection_plane_normals.x_n, resection_plane_normals.y_n ) > 1e-10 && ...
        dot(resection_plane_normals.y_n, resection_plane_normals.z_n) > 1e-10

    disp('Plane normals not perpendicular')

    keyboard

end

% All displacements are defined on the resection plane now based on the
% variable: resection_plane_normals
%% Change position of the cup

% 1) Position on resection surface (superior/inferior, anterior/posterior)
% 2) Offset from resection surface
% 3) Radius (wrt to glenosphere for congruent fit?)
% 4) Depth
% 5) Version/Inclination

% Get mesh data for the hemisphere from Visualisation Object. This will
% then be continuesly updated through "hemisphere" Surface variable

hemi_cup_mesh_data.X = hemisphere_hum.XData;
hemi_cup_mesh_data.Y = hemisphere_hum.YData;
hemi_cup_mesh_data.Z = hemisphere_hum.ZData;

% Cup centre before rototranslation
translated_cup_base = resection_barycentre;
cup_centre_in_humerus = (translated_cup_base + resection_plane_normals.y_n*R)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The vector that will be rotated needs to be a column vector so it is
% transposed here (') then back again after rotations are finished
% Translate cup centre to originate from origin
cup_centre_in_origin = (resection_plane_normals.y_n*R)';

scatter3(cup_centre_in_origin(1), cup_centre_in_origin(2), cup_centre_in_origin(3),'red','filled','o')
scatter3(cup_centre_in_humerus(1), cup_centre_in_humerus(2), cup_centre_in_humerus(3),'cyan','filled','o')

%% Version/Inclination of cup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Axes Rotation Order:
% Z-X-(Y)

% Antero-/Postero- version (about Proximal/Distal axis)
% Rotate Mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Negate (-) hemi_cup_offsets.z_ant_retro_version so that ant/version +ive
hemi_cup_offsets.z_ant_retro_version = - hemi_cup_offsets.z_ant_retro_version;

rotate(hemisphere_hum,...
    resection_plane_normals.z_n,...
    hemi_cup_offsets.z_ant_retro_version,...
    resection_barycentre)

% Rotate cup_centre_in_origin about resection_plane_normals.z_n
R_z = axang2rotm([resection_plane_normals.z_n deg2rad(hemi_cup_offsets.z_ant_retro_version)]);
cup_centre_in_origin = R_z*cup_centre_in_origin;
scatter3(cup_centre_in_origin(1), cup_centre_in_origin(2), cup_centre_in_origin(3),'red','filled','o')

% Need to rotate resection_plane_normals.x_n axis after first rotation
resection_plane_normals.x_n_r1 = R_z*resection_plane_normals.x_n';
resection_plane_normals.x_n_r1 = resection_plane_normals.x_n_r1';

scatter3(resection_plane_normals.x_n_r1(1)*R, resection_plane_normals.x_n_r1(2)*R, resection_plane_normals.x_n_r1(3)*R,'cyan','filled','o')

% Supero-/Infero- inclination (about Anterior/Posterior axis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For this rotation it will be rotated about the new X axis after first Z
% rotation to keep topological meaning for the cup orientation. Can be
% thought of as the orientation of the cup as it sat on the resection plane
% and then rotated about first axes

rotate(hemisphere_hum,...
    resection_plane_normals.x_n_r1,...
    hemi_cup_offsets.x_sup_inf_incl,...  % Rotated X-axis
    resection_barycentre)

% Rotate cup_centre_in_origin about resection_plane_normals.x_n_r1
R_x = axang2rotm([resection_plane_normals.x_n_r1 deg2rad(hemi_cup_offsets.x_sup_inf_incl)]);
cup_centre_in_origin = R_x*cup_centre_in_origin;
scatter3(cup_centre_in_origin(1), cup_centre_in_origin(2), cup_centre_in_origin(3),'red','filled','o')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The vector that was rotated transposed back to row
cup_centre_in_origin = cup_centre_in_origin';

% After rotation translate cup centre in origin back to resection plane
% origin
cup_centre_in_humerus = cup_centre_in_origin + translated_cup_base;
scatter3(cup_centre_in_humerus(1), cup_centre_in_humerus(2), cup_centre_in_humerus(3),'cyan','filled','o')

%% Position on resection surface (anterior/posterior, base offset, superior/inferior)

% % % % Translation offsets in meters (m)
% % % hemi_cup_offsets.x_ant_post   = 0; % X-normal
% % % hemi_cup_offsets.y_base_off   = 0.01; % Y-normal
% % % hemi_cup_offsets.z_prox_dist  = 0.005; % Z-normal

% Anterior / Posterior offsets
hemisphere_hum.XData = hemisphere_hum.XData + resection_plane_normals.x_n(1)*hemi_cup_offsets.x_ant_post;
hemisphere_hum.YData = hemisphere_hum.YData + resection_plane_normals.x_n(2)*hemi_cup_offsets.x_ant_post;
hemisphere_hum.ZData = hemisphere_hum.ZData + resection_plane_normals.x_n(3)*hemi_cup_offsets.x_ant_post;

translated_cup_base = translated_cup_base + resection_plane_normals.x_n*hemi_cup_offsets.x_ant_post;

% Base offset
hemisphere_hum.XData = hemisphere_hum.XData + resection_plane_normals.y_n(1)*hemi_cup_offsets.y_base_off;
hemisphere_hum.YData = hemisphere_hum.YData + resection_plane_normals.y_n(2)*hemi_cup_offsets.y_base_off;
hemisphere_hum.ZData = hemisphere_hum.ZData + resection_plane_normals.y_n(3)*hemi_cup_offsets.y_base_off;

translated_cup_base = translated_cup_base + resection_plane_normals.y_n*hemi_cup_offsets.y_base_off;

% Proximal / Distal offsets
hemisphere_hum.XData = hemisphere_hum.XData + resection_plane_normals.z_n(1)*hemi_cup_offsets.z_prox_dist;
hemisphere_hum.YData = hemisphere_hum.YData + resection_plane_normals.z_n(2)*hemi_cup_offsets.z_prox_dist;
hemisphere_hum.ZData = hemisphere_hum.ZData + resection_plane_normals.z_n(3)*hemi_cup_offsets.z_prox_dist;

translated_cup_base = translated_cup_base + resection_plane_normals.z_n*hemi_cup_offsets.z_prox_dist;


% Calculate position of cup centre aftre rototranslation
% Translate it out with updated cup_base translation offsets
cup_centre_in_humerus = cup_centre_in_origin + translated_cup_base;

% Check that the distance between translated_cup_case and
% cup_centre_in_humerus ewuals the radius
rototran_cup_base_to_centre.vector = translated_cup_base - cup_centre_in_humerus;
rototran_cup_base_to_centre.magnitude = norm(rototran_cup_base_to_centre.vector);

if rototran_cup_base_to_centre.magnitude - R > 1e-4
    disp('Rototranslated cup centre not correct (cup_centre_in_hum - R > 1e-4)')
    keyboard
end

% Rototranslated cup centre and base
scatter3(cup_centre_in_humerus(1), cup_centre_in_humerus(2), cup_centre_in_humerus(3),'green','filled','o')
scatter3(translated_cup_base(1), translated_cup_base(2), translated_cup_base(3),'magenta','filled','o')

% % % keyboard

humerus = struct('resection_plane_normals', resection_plane_normals,...     % x-y-z normal vectors and end-points
    'hemi_cup_mesh_data', hemi_cup_mesh_data,...                            % DEFAULT data for glenoid hemisphere on glenoid plane before any rototranslation
    'hemi_cup_offsets', hemi_cup_offsets,...                                % glenoid hemisphere rototranslation offsets
    'resection_barycentre', resection_barycentre,...                        % centre of DEFAULT hemisphere (CoR of joint effectivly - needs the hemi_gle_offsets translation only to be applied along normals)
    'hemisphere_hum', hemisphere_hum,...                                    % rototranslated cup hemisphere
    'resection_plane', resection_plane,...                                  % resection plane parameters
    'plane_mesh_data', plane_mesh_data,...
    'R', R,...                                                              % hemisphere radius (m)
    'cup_centre_in_humerus', cup_centre_in_humerus,...
    'translated_cup_base', translated_cup_base,...                      % base of humeral cup (should be placed d = R from CoR_glen)
    'stl_hum', stl_hum);

stlwrite_user(['..\..\OpenSim\In\Geometry\cup_' rhash '.stl'],...
    hemisphere_hum.XData,...
    hemisphere_hum.YData,...
    hemisphere_hum.ZData,...
    'mode','ascii',...
    'triangulation','f');

end
