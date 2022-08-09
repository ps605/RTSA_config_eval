function new_model_file = optimDeltViaPoint(model_file)
%optimDeltViaPoint Defines DELT1, DELT2 and DELT3 via poit locations
%   Detailed explanation goes here

%% Set-up
import org.opensim.modeling.*

osim_model = Model(model_file);

% Ackland et al (2010) RTSA MA data

data_RTSA.angles =  [2.5, 30, 60, 90, 120];
data_RTSA.DELT1 =   [15.6, 25.2,32.5, 35.8, 33.3];
data_RTSA.DELT2 =   [30.2, 33.9, 42.2, 46.2, 39.8];
data_RTSA.DELT3 =   [1.3, 3.5, 7.3, 11.4, 14.1];

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

% Some check for it being via point HERE

% Get via point and downcast to get/set location
delt1_via1 = delt1_PPS.get(0);
delt1_via1_downCast = PathPoint.safeDownCast(delt1_via1);

% Get via point location as vector for handling- not Vec3
delt1_via1_loc_i = [delt1_via1_downCast.get_location().get(0),...
    delt1_via1_downCast.get_location().get(1),...
    delt1_via1_downCast.get_location().get(2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELT2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get GeometryPath to calculate MomentArm later
delt2_GP = delt2.getGeometryPath();
delt2_PPS = delt2_GP.getPathPointSet();
delt2_PPS_size = delt2_PPS.getSize();

% Some check for it being via point HERE

% Get via point and downcast to get/set location
delt2_via1 = delt2_PPS.get(0);
delt2_via1_downCast = PathPoint.safeDownCast(delt2_via1);

% Get via point location as vector for handling- not Vec3
delt2_via1_loc_i = [delt2_via1_downCast.get_location().get(0),...
    delt2_via1_downCast.get_location().get(1),...
    delt2_via1_downCast.get_location().get(2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELT3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get GeometryPath to calculate MomentArm later
delt3_GP = delt3.getGeometryPath();
delt3_PPS = delt3_GP.getPathPointSet();
delt3_PPS_size = delt3_PPS.getSize();

% Some check for it being via point HERE

% Get via point and downcast to get/set location
delt3_via1 = delt3_PPS.get(0);
delt3_via1_downCast = PathPoint.safeDownCast(delt3_via1);

% Get via point location as vector for handling- not Vec3
delt3_via1_loc_i = [delt3_via1_downCast.get_location().get(0),...
    delt3_via1_downCast.get_location().get(1),...
    delt3_via1_downCast.get_location().get(2)];

%% Optimisation - MA calculation after via point and joint modulation

% Call to fmincon optim somewhere here
for i_angle = 1:length(data_RTSA.angles)
    
    % Set model to have shoulder_elv value analysed
    osim_model.updCoordinateSet().get('shoulder_elv').setValue(init_state, data_RTSA.angles(i_angle))
    osim_model.realizePosition(init_state)

    % Compute moment arms


end

end