function new_model_file = optimDeltViaPoint(model_file)
%optimDeltViaPoint Defines DELT1, DELT2 and DELT3 via poit locations
%   Detailed explanation goes here

%% Set-up
import org.opensim.modeling.*

osim_model = Model(model_file);
init_state = osim_model.initSystem();

% Ackland et al (2010) RTSA MA data

data_RTSA.angles        = [2.5, 30, 60, 90, 120];
data_RTSA.DELT1_mean    = [15.6, 25.2,32.5, 35.8, 33.3]*0.001;
data_RTSA.DELT1_sd      = [2.3, 2.5, 1.8, 3.4, 2.7]*0.001;
data_RTSA.DELT2_mean    = [30.2, 33.9, 42.2, 46.2, 39.8]*0.001;
data_RTSA.DELT2_sd      = [6.6, 5.5, 5.5, 4.2, 5.8]*0.001;
data_RTSA.DELT3_mean    = [1.3, 3.5, 7.3, 11.4, 14.1]*0.001;
data_RTSA.DELT3_sd      = [1.4, 1.0, 1.7, 2.5, 3.6]*0.001;

%% Handle model

% Get muscle handles
delt1 = osim_model.getMuscles.get('DELT1');
delt2 = osim_model.getMuscles.get('DELT2');
delt3 = osim_model.getMuscles.get('DELT3');

% Get coordinate handle
shoulder_elv = osim_model.getCoordinateSet().get('shoulder_elv');
coord = 'shoulder_elv';

% Get PathPointSet(s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELT1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get GeometryPath to calculate MomentArm later
delt1_GP = delt1.getGeometryPath();
delt1_PPS = delt1_GP.getPathPointSet();
delt1_PPS_size = delt1_PPS.getSize();

% Loop through PathPointSet to get via point
for i_point = 0:delt1_PPS_size-1

    point = delt1_PPS.get(i_point);
    point_name = char(point.getName());

    if strcmp(point_name(end-2:end), 'via')

        % Get via point and downcast to get/set location
        delt1_via_downCast = PathPoint.safeDownCast(point);

        % Get via point location as vector for handling- not Vec3
        delt1_via_loc = [delt1_via_downCast.get_location().get(0),...
            delt1_via_downCast.get_location().get(1),...
            delt1_via_downCast.get_location().get(2)];

    else
        continue
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELT2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get GeometryPath to calculate MomentArm later
delt2_GP = delt2.getGeometryPath();
delt2_PPS = delt2_GP.getPathPointSet();
delt2_PPS_size = delt2_PPS.getSize();

% Some check for it being via point HERE
% Loop through PathPointSet to get via point
for i_point = 0:delt2_PPS_size-1

    point = delt2_PPS.get(i_point);
    point_name = char(point.getName());

    if strcmp(point_name(end-2:end), 'via')

        % Get via point and downcast to get/set location
        delt2_via_downCast = PathPoint.safeDownCast(point);

        % Get via point location as vector for handling- not Vec3
        delt2_via_loc = [delt2_via_downCast.get_location().get(0),...
            delt2_via_downCast.get_location().get(1),...
            delt2_via_downCast.get_location().get(2)];

    else
        continue
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELT3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get GeometryPath to calculate MomentArm later
delt3_GP = delt3.getGeometryPath();
delt3_PPS = delt3_GP.getPathPointSet();
delt3_PPS_size = delt3_PPS.getSize();

for i_point = 0:delt3_PPS_size-1

    point = delt3_PPS.get(i_point);
    point_name = char(point.getName());

    if strcmp(point_name(end-2:end), 'via')

        % Get via point and downcast to get/set location
        delt3_via_downCast = PathPoint.safeDownCast(point);

        % Get via point location as vector for handling- not Vec3
        delt3_via_loc = [delt3_via_downCast.get_location().get(0),...
            delt3_via_downCast.get_location().get(1),...
            delt3_via_downCast.get_location().get(2)];

    else
        continue
    end
end

%% Optimisation - MA calculation after via point and joint modulation
% Search radius around init location
radius = 0.05;
p_sim_0 = delt1_via_loc;

ub = [0.5, 0.5, 0.5];
lb = [-0.5, -0.5, -0.5];

figure(1);
scatter3(p_sim_0(1), p_sim_0(2), p_sim_0(3), 'o', 'filled','magenta')
hold on

% Inequality constraint to keep new via point location with radius of x of
% original via point
f_con = @(p_sim)sphere_func_con(p_sim, p_sim_0, radius);
% Cost function to minimise moment arm differances between simulated and
% calculated conditions 
J = @(p_sim)J_momentArmDist(p_sim, data_RTSA, osim_model, 'DELT1', delt1_via_downCast);

options = optimset('MaxIter', 100, 'TolFun', 1e-6, 'Display', 'iter');

% Run fmincon
[p_sim, fval] = fmincon(J,...
    p_sim_0,...
    [],...
    [],...
    [],...
    [],...
    lb,...
    ub,...
    f_con,...
    options);

figure(1);
scatter3(p_sim_0(1), p_sim_0(2), p_sim_0(3), 'o', 'filled','magenta')
scatter3(p_sim(1), p_sim(2), p_sim(3), 'o', 'filled','red')
xlim([p_sim_0(1)-radius,p_sim_0(1)+radius])
ylim([p_sim_0(2)-radius,p_sim_0(2)+radius])
zlim([p_sim_0(3)-radius,p_sim_0(3)+radius])

figure(2);
scatter3(p_sim_0(1), p_sim_0(2), p_sim_0(3), 'o', 'filled','magenta')
hold on
scatter3(p_sim(1), p_sim(2), p_sim(3), 'o', 'filled','cyan')
% % % Set model to have shoulder_elv value analysed
% % osim_model.updCoordinateSet().get('shoulder_elv').setValue(init_state, deg2rad(data_RTSA.angles(2)));
% % osim_model.realizePosition(init_state);
% % 
% % % %     new_state = osim_model.initSystem()
% % % %     shoulder_elv.getValue(init_state)
% % % Compute moment arms
% % moment_arms.delt1 = delt1_GP.computeMomentArm(init_state, shoulder_elv);
% % moment_arms.delt2 = delt2_GP.computeMomentArm(init_state, shoulder_elv);
% % moment_arms.delt3 = delt3_GP.computeMomentArm(init_state, shoulder_elv);
% % moment_arms


end