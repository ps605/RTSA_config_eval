function J = J_momentArmDist(p_sim, data_abd_RTSA, data_flx_RTSA, osim_model, muscle_name, delt_via_downCast, flag_via_delete )

import org.opensim.modeling.*

init_state = osim_model.initSystem();

% Weights for joint positions
w1      = 10;
w2      = 10;
w3      = 10;
w4      = 10;
w5      = 10;
w6      = 0;
w7      = 0;
w8      = 0;
w9      = 0;
w10     = 0;

% Get coordinate handle
shoulder_elv = osim_model.getCoordinateSet().get('shoulder_elv');
elv_angle    = osim_model.getCoordinateSet().get('elv_angle');

% Get muscle handles
deltx = osim_model.getMuscles.get(muscle_name);

% Get GeometryPaths
if strcmp(muscle_name, 'DELT1')
    delt_GP = deltx.getGeometryPath();
    exp_abd_MA_mean = data_abd_RTSA.DELT1_mean;
    exp_flx_MA_mean = data_flx_RTSA.DELT1_mean;
elseif strcmp(muscle_name, 'DELT2') 
    delt_GP = deltx.getGeometryPath();
    exp_abd_MA_mean = data_abd_RTSA.DELT2_mean;
    exp_flx_MA_mean = data_flx_RTSA.DELT2_mean;
elseif strcmp(muscle_name, 'DELT3')
    delt_GP = deltx.getGeometryPath();
    exp_abd_MA_mean = data_abd_RTSA.DELT3_mean;
    exp_flx_MA_mean = data_flx_RTSA.DELT3_mean;
    
end

if flag_via_delete == false
    % Set new location for DeltX via point
    delt_via_downCast.set_location(Vec3(p_sim(1), p_sim(2), p_sim(3)))
    % Call ::Model.finalizeConnections() to....
    osim_model.finalizeConnections();
    % Re initialise the system to allow computing moment arm
    new_state = osim_model.initSystem;
else
    new_state = osim_model.initSystem;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 1 - 2.5 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
osim_model.updCoordinateSet().get('shoulder_elv').setValue(new_state, deg2rad(data_abd_RTSA.angles(1)), true);
osim_model.realizePosition(new_state);


delt_abd_MA.pos1 = delt_GP.computeMomentArm(new_state, shoulder_elv);
delt_abd_MA.diff1 = delt_abd_MA.pos1 - exp_abd_MA_mean(1);

delt_flx_MA.pos1 = delt_GP.computeMomentArm(new_state, elv_angle);
delt_flx_MA.diff1 = delt_flx_MA.pos1 - exp_flx_MA_mean(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 2 - 30 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
osim_model.updCoordinateSet().get('shoulder_elv').setValue(new_state, deg2rad(data_abd_RTSA.angles(2)), true);
osim_model.realizePosition(new_state);

delt_abd_MA.pos2 = delt_GP.computeMomentArm(new_state, shoulder_elv);
delt_abd_MA.diff2 = delt_abd_MA.pos2 - exp_abd_MA_mean(2);

delt_flx_MA.pos2 = delt_GP.computeMomentArm(new_state, elv_angle);
delt_flx_MA.diff2 = delt_flx_MA.pos2 - exp_flx_MA_mean(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 3 - 60 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
osim_model.updCoordinateSet().get('shoulder_elv').setValue(new_state, deg2rad(data_abd_RTSA.angles(3)), true);
osim_model.realizePosition(new_state);

delt_abd_MA.pos3 = delt_GP.computeMomentArm(new_state, shoulder_elv);
delt_abd_MA.diff3 = delt_abd_MA.pos3 - exp_abd_MA_mean(3);

delt_flx_MA.pos3 = delt_GP.computeMomentArm(new_state, elv_angle);
delt_flx_MA.diff3 = delt_flx_MA.pos3 - exp_flx_MA_mean(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 4 - 90 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
osim_model.updCoordinateSet().get('shoulder_elv').setValue(new_state, deg2rad(data_abd_RTSA.angles(4)), true);
osim_model.realizePosition(new_state);

delt_abd_MA.pos4 = delt_GP.computeMomentArm(new_state, shoulder_elv);
delt_abd_MA.diff4 = delt_abd_MA.pos4 - exp_abd_MA_mean(4);

delt_flx_MA.pos4 = delt_GP.computeMomentArm(new_state, elv_angle);
delt_flx_MA.diff4 = delt_flx_MA.pos4 - exp_flx_MA_mean(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 5 - 120 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
osim_model.updCoordinateSet().get('shoulder_elv').setValue(new_state, deg2rad(data_abd_RTSA.angles(5)), true);
osim_model.realizePosition(new_state);

delt_abd_MA.pos5 = delt_GP.computeMomentArm(new_state, shoulder_elv);
delt_abd_MA.diff5 = delt_abd_MA.pos5 - exp_abd_MA_mean(5);

delt_flx_MA.pos5 = delt_GP.computeMomentArm(new_state, elv_angle);
delt_flx_MA.diff5 = delt_flx_MA.pos5 - exp_flx_MA_mean(5);

% Assemble model
osim_model.assemble(new_state)
%% Calculate sum of moment arm difference 
J = (abs(w1*delt_abd_MA.diff1) + abs(w2*delt_abd_MA.diff2) + abs(w3*delt_abd_MA.diff3) + abs(w4*delt_abd_MA.diff4) + abs(w5*delt_abd_MA.diff5) +...
    abs(w6*delt_flx_MA.diff1) + abs(w7*delt_flx_MA.diff2) + abs(w8*delt_flx_MA.diff3) + abs(w9*delt_flx_MA.diff4) + abs(w10*delt_flx_MA.diff5))/(w1 + w2 + w3+ w4 + w5 + w6 + w7 + w8 + w9 + w10);
% J = 1000*(sqrt(w1*delt_MA.diff1^2 + w2*delt_MA.diff2^2 + w3*delt_MA.diff3^2 + w4*delt_MA.diff4^2 + w5*delt_MA.diff5^2))/(w1 + w2 + w3+ w4 + w5);

end