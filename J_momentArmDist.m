function J = J_momentArmDist(p_sim, data_RTSA, osim_model, muscle_name, delt_via_downCast )

import org.opensim.modeling.*

init_state = osim_model.initSystem();

% Weights for joint positions
w1  = 1;
w2  = 0;
w3  = 0;
w4  = 0;
w5  = 0;

% Get coordinate handle
shoulder_elv = osim_model.getCoordinateSet().get('shoulder_elv');

% Get muscle handles
deltx = osim_model.getMuscles.get(muscle_name);

figure(101);
scatter3(p_sim(1), p_sim(2), p_sim(3), 'o', 'filled','cyan')
hold on

% Get GeometryPaths
if strcmp(muscle_name, 'DELT1')
    delt_GP = deltx.getGeometryPath();
    exp_MA_mean = data_RTSA.DELT1_mean;
    exp_MA_sd   = data_RTSA.DELT1_sd;
elseif strcmp(muscle_name, 'DELT2') 
    delt_GP = deltx.getGeometryPath();
    exp_MA_mean = data_RTSA.DELT2_mean;
    exp_MA_sd   = data_RTSA.DELT2_sd;
elseif strcmp(muscle_name, 'DELT3')
    delt_GP = deltx.getGeometryPath();
    exp_MA_mean = data_RTSA.DELT3_mean;
    exp_MA_sd   = data_RTSA.DELT3_sd;
end

% Set new location for DeltX via point
delt_via_downCast.set_location(Vec3(p_sim(1), p_sim(2), p_sim(3)))
% Call ::Model.finalizeConnections() to....
osim_model.finalizeConnections();
% Re initialise the system to allow computing moment arm
new_state = osim_model.initSystem;

%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 1 - 2.5 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
osim_model.updCoordinateSet().get('shoulder_elv').setValue(new_state, deg2rad(data_RTSA.angles(1)), false);
osim_model.realizePosition(new_state);


delt_MA.pos1 = delt_GP.computeMomentArm(new_state, shoulder_elv);
delt_MA.diff1 = delt_MA.pos1 - exp_MA_mean(1);

% % figure(3);
% % scatter(data_RTSA.angles(1), delt_MA.pos1,'filled','o','cyan');
% % hold on

% if delt_MA.pos1 >= exp_MA_mean(1) - exp_MA_sd(1) && delt_MA.pos1 <= exp_MA_mean(1) + exp_MA_sd(1)
%     disp('lol')
%     %keyboard
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 2 - 30 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
osim_model.updCoordinateSet().get('shoulder_elv').setValue(new_state, deg2rad(data_RTSA.angles(2)), false);
osim_model.realizePosition(new_state);

delt_MA.pos2 = delt_GP.computeMomentArm(new_state, shoulder_elv);
delt_MA.diff2 = delt_MA.pos2 - exp_MA_mean(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 3 - 60 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
osim_model.updCoordinateSet().get('shoulder_elv').setValue(new_state, deg2rad(data_RTSA.angles(3)), false);
osim_model.realizePosition(new_state);

delt_MA.pos3 = delt_GP.computeMomentArm(new_state, shoulder_elv);
delt_MA.diff3 = delt_MA.pos3 - exp_MA_mean(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 4 - 90 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
osim_model.updCoordinateSet().get('shoulder_elv').setValue(new_state, deg2rad(data_RTSA.angles(4)), false);
osim_model.realizePosition(new_state);

delt_MA.pos4 = delt_GP.computeMomentArm(new_state, shoulder_elv);
delt_MA.diff4 = delt_MA.pos4 - exp_MA_mean(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 5 - 120 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
osim_model.updCoordinateSet().get('shoulder_elv').setValue(new_state, deg2rad(data_RTSA.angles(5)), false);
osim_model.realizePosition(new_state);

delt_MA.pos5 = delt_GP.computeMomentArm(new_state, shoulder_elv);
delt_MA.diff5 = delt_MA.pos5 - exp_MA_mean(5);

% Assemble model
osim_model.assemble(new_state)
%% Calculate sum of moment arm difference 
J = (abs(w1*delt_MA.diff1) + abs(w2*delt_MA.diff2) + abs(w3*delt_MA.diff3) + abs(w4*delt_MA.diff4) + abs(w5*delt_MA.diff5))/(w1 + w2 + w3+ w4 + w5);
% J = 1000*(sqrt(w1*delt_MA.diff1^2 + w2*delt_MA.diff2^2 + w3*delt_MA.diff3^2 + w4*delt_MA.diff4^2 + w5*delt_MA.diff5^2))/(w1 + w2 + w3+ w4 + w5);

% figure(102);
% title('J')
% scatter(1, J, 'o', 'filled', 'cyan')
% ylim([-0.01 0.05])
% hold on

end