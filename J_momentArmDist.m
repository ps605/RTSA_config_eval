function J = J_momentArmDist(p_sim, data_RTSA, osim_model, muscle_name, delt_via_downCast )

import org.opensim.modeling.*

init_state = osim_model.initSystem();

% Get coordinate handle
shoulder_elv = osim_model.getCoordinateSet().get('shoulder_elv');

% Get muscle handles
deltx = osim_model.getMuscles.get(muscle_name);

figure(1);
scatter3(p_sim(1), p_sim(2), p_sim(3), 'o', 'filled','cyan')
hold on

% Get GeometryPaths
if strcmp(muscle_name, 'DELT1')
    delt_GP = deltx.getGeometryPath();
    exp_MA = data_RTSA.DELT1;
    
elseif strcmp(muscle_name, 'DELT2') 
    delt_GP = deltx.getGeometryPath();
    exp_MA = data_RTSA.DELT2;
elseif strcmp(muscle_name, 'DELT3')
    delt_GP = deltx.getGeometryPath();
    exp_MA = data_RTSA.DELT3;
end

% Set new location for DeltX via point
delt_via_downCast.set_location(Vec3(p_sim(1), p_sim(2), p_sim(3)))

%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 1 - 2.5 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
osim_model.updCoordinateSet().get('shoulder_elv').setValue(init_state, deg2rad(data_RTSA.angles(1)));
osim_model.realizePosition(init_state);


delt_MA.pos1 = delt_GP.computeMomentArm(init_state, shoulder_elv);
delt_MA.diff1 = delt_MA.pos1 - exp_MA(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 2 - 30 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
osim_model.updCoordinateSet().get('shoulder_elv').setValue(init_state, deg2rad(data_RTSA.angles(2)));
osim_model.realizePosition(init_state);

delt_MA.pos2 = delt_GP.computeMomentArm(init_state, shoulder_elv);
delt_MA.diff2 = delt_MA.pos2 - exp_MA(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 3 - 60 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
osim_model.updCoordinateSet().get('shoulder_elv').setValue(init_state, deg2rad(data_RTSA.angles(3)));
osim_model.realizePosition(init_state);

delt_MA.pos3 = delt_GP.computeMomentArm(init_state, shoulder_elv);
delt_MA.diff3 = delt_MA.pos3 - exp_MA(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 4 - 90 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
osim_model.updCoordinateSet().get('shoulder_elv').setValue(init_state, deg2rad(data_RTSA.angles(4)));
osim_model.realizePosition(init_state);

delt_MA.pos4 = delt_GP.computeMomentArm(init_state, shoulder_elv);
delt_MA.diff4 = delt_MA.pos4 - exp_MA(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 5 - 120 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
osim_model.updCoordinateSet().get('shoulder_elv').setValue(init_state, deg2rad(data_RTSA.angles(5)));
osim_model.realizePosition(init_state);

delt_MA.pos5 = delt_GP.computeMomentArm(init_state, shoulder_elv);
delt_MA.diff5 = delt_MA.pos5 - exp_MA(5);

%% Calculate sum of moment arm difference 

J = delt_MA.diff1^2 + delt_MA.diff2^2 + delt_MA.diff3^2 + delt_MA.diff4^2 + delt_MA.diff5^2;

end