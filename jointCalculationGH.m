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
