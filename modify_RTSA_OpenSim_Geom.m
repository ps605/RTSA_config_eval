% Parametric modifications of RTSA implant geometry configurations (humeral cup and
% glenosphere) and shoulder anatomy (humerus and scapula) to approximate
% "Virtual Surgery" in RTSA. Main output of this program is to adjust the
% relative position of the humerus with respect to the scapula after the
% RTSA implant components with selected configurations have been positioned
% on the respective anatomies. The configurations include congruent cup and
% hemisphere radii (R_cup = R_hemisphere), position (antero/posterior,
% supero/inferial and base offset) and orientation (antero/postero,
% supero/infero version) of the components. When geometric modification
% have been completed on each anatomy and the implant component positions
% are defined the parametic approximation of the components are exported as
% .stl and then the two anatomies and implant componets are registered in
% the global coordinate system. The newly created shoulder joint CoR in the
% scapula (parent body) and humerus (child body) are exported and redefined
% in OpenSim.
%
% Critical functions:
% axang2rotm, vrrotvec, rotate, fmincon
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% TO ADD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Reaming depth on glenoid (translate glenoid plane along -Z norm of
%       glenoid plane
% 2) Interactivly chose surgical configuration parameters and updating
%       visualiser and data to pass onto next step?
%       (https://www.mathworks.com/help/control/ug/build-app-with-interactive-plot-updates.html)
% 3) Dynamic reading of .stl files (ALMOST - 20220630)
% 4) Where the hemisphere and cup .stl are saved to 
% 5) Link final GHJ outputs to OpenSim to create new model (COMPLETE
%   20220629) - now to run simulation
% 
%
% Pavlos Silvestros, PhD - University of Victoria, CAN, June 2022

tic
close all;
clear;
clc;

%% Define Parameters for hemisphere/cup gemetry and offsets
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hemisphere radius %%%%%%%%%%%%%%%%%%%%%%%%%%%%
diameter = 0.036;

R = diameter/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Glenosphere offsets %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rotation offsets in degrees

% Anteroversion: +ive; Retroversion: -ive
hemi_gle_offsets.y_ant_retro_version    = 0;
% Inferior inclination: - ive; Superior inclination: +ive
hemi_gle_offsets.x_sup_inf_incl         = -10;

% Translation offsets in meters (m)
hemi_gle_offsets.x_ant_post   = 0;          % X-normal
hemi_gle_offsets.y_prox_dist  = -0.006;     % Y-normal
hemi_gle_offsets.z_base_off   = 0.006;      % Z-normal

%%%%%%%%%%%%%%%%%%%%%%%%%%% Humeral cup offsets %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rotation offsets in degrees

% Anteroversion: +ive; Retroversion: -ive
hemi_cup_offsets.z_ant_retro_version   = 0;
% Inferior inclination: - ive; Superior inclination: +ive
hemi_cup_offsets.x_sup_inf_incl        = 12.5;

% Translation offsets in meters (m)
hemi_cup_offsets.x_ant_post   = 0;      % X-normal
hemi_cup_offsets.y_base_off   = 0.012;  % Y-normal
hemi_cup_offsets.z_prox_dist  = 0;      % Z-normal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Flags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If should plot intermediate plots for checking
flag_checkPlots = false;

% If replacing Muscles with Actuators
flag_useTorque = false;

% If removing Rotator Cuff muscles
flag_keepRC =  false;

% Replace muscle models Millard2012Equilibrium with DeGrootFregly
flag_ReplaceMuscles = true;

% Run Moco after model is defined?
flag_runSim = true;

% Create a random 11-char hash to reference model file X00yyy111zz (~30e12)
rhash = [char(randi([65 90],1,1))...    
    char(randi([48 57],1,2))...         
    char(randi([97 122],1,3))...       
    char(randi([48 57],1,3))...         
    char(randi([97 122],1,2))];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Call functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create internal functions here (one for humerus and one for scapula) that
% plot and do all the positioning then only return necessary values and
% data.

% Define parametric implant on .stl anatomy & extract parameters in global
scapula = glenoidGeom(R, hemi_gle_offsets, rhash);

% Define parametric implant on .stl anatomy & extract parameters in global
humerus = humerusGeom(R, hemi_cup_offsets, rhash);

% Read in defined implant parameters and .stl and calculate GHJ centre
[GHJ_in_parent, GHJ_in_child] = jointCalculationGH(scapula,humerus);

% Define OpenSim model with new GHJ parameters from 'Virtual Surgery'
model_file = adjustOpenSimModelGHJ(GHJ_in_parent,...
    GHJ_in_child,...
    hemi_gle_offsets,...
    hemi_cup_offsets,...
    R,...
    rhash,...
    flag_useTorque,...
    flag_keepRC,...
    flag_ReplaceMuscles);

% Run OpenSim moco for predictive simulation
if flag_runSim == true
    runRTSAsims(model_file, rhash, flag_keepRC)
end

toc

%% Calculate Joint frame positions with OpenSim definitions
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [GHJ_in_parent, GHJ_in_child] = jointCalculationGH(scapula,humerus)

% Everything is positioned wrt global that coincident with Scapula frame
% (fixed) and initial Humeral frame which will be transformed based on GHJ
% positions then all will be plotted

% GHJ in Parent (Scapula)
% This is the position of glenoid hemisphere centre (CoR) from (0, 0 ,0)
GHJ_in_parent = scapula.CoR_glen;

% Offset of the Humeral cup centre in humerus (in inertial frame) to GHJ
cup_centre_to_GHJ = GHJ_in_parent - humerus.cup_centre_in_humerus;
% Apply offset to cup barrycentre - this is now the coincident
% connect_point
cup_base_registered = humerus.translated_cup_base + cup_centre_to_GHJ;

% Apply offset to Humerus origin - (0, 0, 0)
humerus_origin_registered = [0 0 0] + cup_centre_to_GHJ;

% Need to project point on hemisphere surface distance or R from centroid
% to then place the base/centre of cup

% % Select point on hemisphere to conenct with cup - IS THIS NEEDED NOW?
% keyboard
connect_point = cup_base_registered; %cursor_info.Position;

% Check that the norm of the position vector is the same as the hemisphere
% radius
% Calculate position vector from GHJ_in_parent to GHJ_in_child
connect_point_to_CoR_glen.vector = connect_point - scapula.CoR_glen;
% Calculate the magnitude
connect_point_to_CoR_glen.magnitude = norm(connect_point_to_CoR_glen.vector);
% Normalise position to unit vector
connect_point_to_CoR_glen.dir_vector = connect_point_to_CoR_glen.vector/connect_point_to_CoR_glen.magnitude;
if abs(scapula.R - connect_point_to_CoR_glen.magnitude) > 1e-4
    disp('Hemisphere point to CoR normal not equal to radius');
end

% GHJ in Child (Humerus)
% This assumes that the Humerus and attached cup will be translated from
% the global (base) coordinate system as that is where the initial
% calculations and exported data were done. Now they will be translated
% down (by original cup centre to GHJ_in_parent [i.e. centre to centre] so
% that the cup base conects to the 'connect_point')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Back off from the connect_point(registered translated_cup_base point from centre_to_centre
% translation) by the same offset as the translated barycentre in original
% humerus centre. This is then the position of the Humerus body origin in
% the Scapula after which the GHJ_in_child can be calculated.

% As cup is connected ridgidly to humerus the same translation is true for
% the humerus reference frame from origin (0, 0, 0)

% Position of Humerus reference frame in Scapula is equal to centre to
% centre offset (i.e. humerus_centre_registered)

GHJ_in_child = GHJ_in_parent - humerus_origin_registered;

line([GHJ_in_parent(1) humerus_origin_registered(1)],...
    [GHJ_in_parent(2) humerus_origin_registered(2)],...
    [GHJ_in_parent(3) humerus_origin_registered(3)], ...
    'LineWidth',4,'Color','magenta');

%% Plot Scapula stl (This will be the base for the humerus also)

figure(20);

% Plot global coordinate system
x_hat=[0.1 0 0];
y_hat=[0 0.1 0];
z_hat=[0 0 0.1];

line([x_hat(1) 0],[x_hat(2) 0],[x_hat(3) 0], 'LineWidth',4,'Color','r'); % X - Red
line([y_hat(1) 0],[y_hat(2) 0],[y_hat(3) 0], 'LineWidth',4,'Color','y'); % Y - Yellow
line([z_hat(1) 0],[z_hat(2) 0],[z_hat(3) 0], 'LineWidth',4,'Color','g'); % Z - Green

% Plot .stl
patch(scapula.stl_scap.x, scapula.stl_scap.y, scapula.stl_scap.z,'b',...
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

%% Plot Glenoid plane

% Plot Plane
surf(scapula.plane_mesh_data.x_plane, scapula.plane_mesh_data.y_plane, scapula.plane_mesh_data.z_plane,...
    'FaceColor','g',...
    'FaceAlpha', 0.5,...
    'EdgeAlpha', 0.25)

%% Plot Glenoid hemisphere, CoR and coordinate system

% Default barycentre
scatter3(scapula.glenoid_barycentre(1), scapula.glenoid_barycentre(2), scapula.glenoid_barycentre(3), 'magenta', 'filled', 'o');

% The hemisphere with the CoR
surf(scapula.hemisphere_gle.XData, scapula.hemisphere_gle.YData, scapula.hemisphere_gle.ZData,...
    'FaceColor',[ 1 1 0],...
    'FaceAlpha', 0.75,...
    'EdgeColor', [0 0 0],...
    'EdgeAlpha', 0.1);

% The CoR from the moved hemisphere
scatter3(scapula.CoR_glen(1), scapula.CoR_glen(2), scapula.CoR_glen(3), 'magenta', 'filled', 'o')

% The plane normals for the glenoid

% X - Anterior/Posterior
scatter3(scapula.glenoid_plane_normals.x_p(1), scapula.glenoid_plane_normals.x_p(2), scapula.glenoid_plane_normals.x_p(3),'red','filled','o')
line([scapula.glenoid_barycentre(1) scapula.glenoid_plane_normals.x_p(1)],...
    [scapula.glenoid_barycentre(2) scapula.glenoid_plane_normals.x_p(2)],...
    [scapula.glenoid_barycentre(3) scapula.glenoid_plane_normals.x_p(3)], ...
    'LineWidth',4,'Color','red');

% Y - Proximal/Distal
scatter3(scapula.glenoid_plane_normals.y_p(1), scapula.glenoid_plane_normals.y_p(2), scapula.glenoid_plane_normals.y_p(3),'yellow','filled','o')
line([scapula.glenoid_barycentre(1) scapula.glenoid_plane_normals.y_p(1)],...
    [scapula.glenoid_barycentre(2) scapula.glenoid_plane_normals.y_p(2)],...
    [scapula.glenoid_barycentre(3) scapula.glenoid_plane_normals.y_p(3)], ...
    'LineWidth',4,'Color','yellow');

% Z - Out of plane
scatter3(scapula.glenoid_plane_normals.z_p(1), scapula.glenoid_plane_normals.z_p(2), scapula.glenoid_plane_normals.z_p(3),'green','filled','o')
line([scapula.glenoid_barycentre(1) scapula.glenoid_plane_normals.z_p(1)],...
    [scapula.glenoid_barycentre(2) scapula.glenoid_plane_normals.z_p(2)],...
    [scapula.glenoid_barycentre(3) scapula.glenoid_plane_normals.z_p(3)], ...
    'LineWidth',4,'Color','green');



figure(20)
% Plot connecting lines (translation checks)
% Line connecting cup base when Humerus at (0, 0, 0) to connect_point (i.e.
% the translation needed to translate the Humerus at GHJ)
line([humerus.translated_cup_base(1) connect_point(1)],...
    [humerus.translated_cup_base(2) connect_point(2)],...
    [humerus.translated_cup_base(3) connect_point(3)], ...
    'LineWidth',4,'Color','magenta');
% Visualise cup base in Humerus (origin)
scatter3(humerus.translated_cup_base(1),humerus.translated_cup_base(2),humerus.translated_cup_base(3),'magenta','filled','o')
% Visualise to check (connect_point = cup_base_registered);
scatter3(connect_point(1),connect_point(2),connect_point(3),'black','filled','o')

% Visualise registered Humeral cup base on hemisphere (from centre to
% centre offset)
scatter3(cup_base_registered(1),cup_base_registered(2),cup_base_registered(3),'green','filled','o')

% Visualise (Humerus cup centre in origin)
scatter3(humerus.cup_centre_in_humerus(1),humerus.cup_centre_in_humerus(2),humerus.cup_centre_in_humerus(3),'magenta','filled','o')


% Coincident connect_point (hemi and cup) to CoR
line([connect_point(1) GHJ_in_parent(1)],...
    [connect_point(2) GHJ_in_parent(2)],...
    [connect_point(3) GHJ_in_parent(3)], ...
    'LineWidth',4,'Color','cyan');
% Visualise the CoR (hemisphere centre)
scatter3(GHJ_in_parent(1), GHJ_in_parent(2), GHJ_in_parent(3),'black','filled','o')

% Cup centre to GHJ_in_parent
line([humerus.cup_centre_in_humerus(1) GHJ_in_parent(1)],...
    [humerus.cup_centre_in_humerus(2) GHJ_in_parent(2)],...
    [humerus.cup_centre_in_humerus(3) GHJ_in_parent(3)], ...
    'LineWidth',4,'Color','magenta');

% Origin (0, 0, 0)
scatter3(0, 0, 0,'black','filled','o');

% Visualise vector between the registered Humeral origin to the CoR
line([humerus_origin_registered(1) GHJ_in_parent(1)],...
    [humerus_origin_registered(2) GHJ_in_parent(2)],...
    [humerus_origin_registered(3) GHJ_in_parent(3)], ...
    'LineWidth',4,'Color','green');
% Visualise the CoR (hemisphere centre)
scatter3(humerus_origin_registered(1), humerus_origin_registered(2), humerus_origin_registered(3),'black','filled','o')

% Plot original Humeral cup in Origin
surf(humerus.hemisphere_hum.XData,...
    humerus.hemisphere_hum.YData,...
    humerus.hemisphere_hum.ZData,...
    'FaceColor',[ 0 0 1],...
    'FaceAlpha', 0.85,...
    'EdgeColor', [0 0 0],...
    'EdgeAlpha', 0.1);


% Plot registered cup. Copy XYZ data from object and translate the offset
% between cup centre in humerus (original) to GHJ (centre to centre)

cup_hum_trans.XData = humerus.hemisphere_hum.XData + cup_centre_to_GHJ(1);
cup_hum_trans.YData = humerus.hemisphere_hum.YData + cup_centre_to_GHJ(2);
cup_hum_trans.ZData = humerus.hemisphere_hum.ZData + cup_centre_to_GHJ(3);

humerus.hemisphere_hum_reg = surf(cup_hum_trans.XData,...
    cup_hum_trans.YData,...
    cup_hum_trans.ZData,...
    'FaceColor',[ 1 0 0],...
    'FaceAlpha', 0.85,...
    'EdgeColor', [0 0 0],...
    'EdgeAlpha', 0.1);

% Plot registered humerus .stl
patch(humerus.stl_hum.x + cup_centre_to_GHJ(1),...
    humerus.stl_hum.y + cup_centre_to_GHJ(2),...
    humerus.stl_hum.z + cup_centre_to_GHJ(3),...
    'b',...
    'FaceColor', [0.8 0.8 0.8],...
    'FaceAlpha', 0.75,...
    'EdgeColor', [0.65 0.65 0.65],...
    'EdgeAlpha', 0.5);

%% Plot clean final figure with registered geometries

figure(101)
title('Registered Scapular and Humeral geometries with parametric implant configuration after "Virtual Surgery"')

% Plot global coordinate system
x_hat=[0.1 0 0];
y_hat=[0 0.1 0];
z_hat=[0 0 0.1];

line([x_hat(1) 0],[x_hat(2) 0],[x_hat(3) 0], 'LineWidth',4,'Color','r'); % X - Red
line([y_hat(1) 0],[y_hat(2) 0],[y_hat(3) 0], 'LineWidth',4,'Color','y'); % Y - Yellow
line([z_hat(1) 0],[z_hat(2) 0],[z_hat(3) 0], 'LineWidth',4,'Color','g'); % Z - Green

% Plot scapula .stl
patch(scapula.stl_scap.x, scapula.stl_scap.y, scapula.stl_scap.z,'b',...
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

% The hemisphere with the CoR
surf(scapula.hemisphere_gle.XData, scapula.hemisphere_gle.YData, scapula.hemisphere_gle.ZData,...
    'FaceColor',[ 1 1 0],...
    'FaceAlpha', 0.75,...
    'EdgeColor', [0 0 0],...
    'EdgeAlpha', 0.1);

% The CoR from the moved hemisphere
scatter3(scapula.CoR_glen(1), scapula.CoR_glen(2), scapula.CoR_glen(3), 'magenta', 'filled', 'o')

% The plane normals for the glenoid

% % % % X - Anterior/Posterior
% % % scatter3(scapula.glenoid_plane_normals.x_p(1), scapula.glenoid_plane_normals.x_p(2), scapula.glenoid_plane_normals.x_p(3),'red','filled','o')
% % % line([scapula.glenoid_barycentre(1) scapula.glenoid_plane_normals.x_p(1)],...
% % %     [scapula.glenoid_barycentre(2) scapula.glenoid_plane_normals.x_p(2)],...
% % %     [scapula.glenoid_barycentre(3) scapula.glenoid_plane_normals.x_p(3)], ...
% % %     'LineWidth',4,'Color','red');
% % %
% % % % Y - Proximal/Distal
% % % scatter3(scapula.glenoid_plane_normals.y_p(1), scapula.glenoid_plane_normals.y_p(2), scapula.glenoid_plane_normals.y_p(3),'yellow','filled','o')
% % % line([scapula.glenoid_barycentre(1) scapula.glenoid_plane_normals.y_p(1)],...
% % %     [scapula.glenoid_barycentre(2) scapula.glenoid_plane_normals.y_p(2)],...
% % %     [scapula.glenoid_barycentre(3) scapula.glenoid_plane_normals.y_p(3)], ...
% % %     'LineWidth',4,'Color','yellow');
% % %
% % % % Z - Out of plane
% % % scatter3(scapula.glenoid_plane_normals.z_p(1), scapula.glenoid_plane_normals.z_p(2), scapula.glenoid_plane_normals.z_p(3),'green','filled','o')
% % % line([scapula.glenoid_barycentre(1) scapula.glenoid_plane_normals.z_p(1)],...
% % %     [scapula.glenoid_barycentre(2) scapula.glenoid_plane_normals.z_p(2)],...
% % %     [scapula.glenoid_barycentre(3) scapula.glenoid_plane_normals.z_p(3)], ...
% % %     'LineWidth',4,'Color','green');

% Plot registered humeral cup
humerus.hemisphere_hum_reg = surf(cup_hum_trans.XData,...
    cup_hum_trans.YData,...
    cup_hum_trans.ZData,...
    'FaceColor',[ 1 0 0],...
    'FaceAlpha', 0.85,...
    'EdgeColor', [0 0 0],...
    'EdgeAlpha', 0.1);

% Plot registered humerus .stl
patch(humerus.stl_hum.x + cup_centre_to_GHJ(1),...
    humerus.stl_hum.y + cup_centre_to_GHJ(2),...
    humerus.stl_hum.z + cup_centre_to_GHJ(3),...
    'b',...
    'FaceColor', [0.8 0.8 0.8],...
    'FaceAlpha', 0.75,...
    'EdgeColor', [0.65 0.65 0.65],...
    'EdgeAlpha', 0.5);
end

%% For the humerus
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function humerus = humerusGeom(R, hemi_cup_offsets,rhash)
%% Set up
% Load in and configure points of Humerurs .stl
[x, y, z] = stlread('..\..\OpenSim\In\Geometry\humerus_resected_manifold_closed.stl');

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

figure (2);
pcshow(resection_pointCloud, 'MarkerSize', 100)

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
    disp('You fucked up')
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

%% To reset run this


    function resetHumerus
        hemisphere_hum.XData = hemi_cup_mesh_data.X;
        hemisphere_hum.YData = hemi_cup_mesh_data.Y;
        hemisphere_hum.ZData = hemi_cup_mesh_data.Z;

        cup_centre_in_origin = (resection_plane_normals.y_n*R)';
    end
end

%% For the scapula
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function scapula = glenoidGeom(R, hemi_gle_offsets, rhash)
%% Set up
% Load in and configure points of Humerurs .stl
[x, y, z] = stlread('..\..\OpenSim\In\Geometry\scapula.stl');

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

% Pause to select points on .stl where "resection" was made and load in from
% plot
% % % keyboard;

% Get from base workspace because the object exports data tip there...
% glenoid_points = evalin('base', 'cursor_info');

p = load('glenoid_points_checked.mat');
%% Fit plane to glenoid rim points

% Glenoid rim points X-Y-Z
glenoid_points = vertcat(p.glenoid_points(:).Position);
% Calculate the mean of the points
glenoid_barycentre = mean(glenoid_points);

% Get points to define coordinate system - Y-axis
[~, min_gl_p] = min(glenoid_points(:,2));
glenoid_inferior = glenoid_points(min_gl_p,:);

% Calculate vector between barrycentre and most inferior point
bc_to_inferior_p =  glenoid_inferior - glenoid_barycentre;


%% continue
% Invert so it point superiorly and normalise
axial_normal = -(bc_to_inferior_p/norm(bc_to_inferior_p));

scatter3(glenoid_points(:,1), glenoid_points(:,2), glenoid_points(:,3), 'filled', 'o', 'cyan');
scatter3(glenoid_points(min_gl_p,1), glenoid_points(min_gl_p,2), glenoid_points(min_gl_p,3), 'filled', 'o', 'magenta');

% Generate PointCloud of slected points and barycentre
glenoid_pointCloud = pointCloud(vertcat(glenoid_points, glenoid_barycentre));

figure (2);
pcshow(glenoid_pointCloud, 'MarkerSize', 100)

% Fit plane to the points
[glenoid_plane,~,~, ~] = pcfitplane(glenoid_pointCloud, 0.0001);
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
[glenoid_inferior_on_plane,fval] = fmincon(J_inferior,...
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
[glenoid_barycentre,fval] = fmincon(J_barycentre,...
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

sca_bary = [x_sca_mean y_sca_mean z_sca_mean];

scatter3(sca_bary(1), sca_bary(2), sca_bary(3),'cyan', 'filled','o')

% Pass .stl x-y-z to structure to then export
stl_scap = struct('x', x,...
    'y', y,...
    'z', z);

% What is the angle between the norm and the vector from the norm origine to glenoid vertex barycentre

% Vector from glenoid barycentre to humeral verteces barrycentre
gle_bary_vec = sca_bary - glenoid_barycentre;

% Connect with line to visualise vector
line([glenoid_barycentre(1) sca_bary(1)],...
    [glenoid_barycentre(2) sca_bary(2)],...
    [glenoid_barycentre(3) sca_bary(3)], ...
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

%% To reset run this
    function resetGlenoid
        hemisphere_gle.XData = hemi_gle_mesh_data.X;
        hemisphere_gle.YData = hemi_gle_mesh_data.Y;
        hemisphere_gle.ZData = hemi_gle_mesh_data.Z;

        CoR_glen = glenoid_barycentre;
    end



end
