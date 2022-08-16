function plotViaOptResults(muscle_name, p_sim_0, p_sim_opt, delt_via_downCast, delt_GP, shoulder_elv, osim_model, data_RTSA, model_MA_init, J_opt, radius)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


import org.opensim.modeling.*

init_state = osim_model.initSystem();

% Select which DELTx
if strcmp(muscle_name, 'DELT1')
    exp_MA_mean = data_RTSA.DELT1_mean;
    exp_MA_sd   = data_RTSA.DELT1_sd;
elseif strcmp(muscle_name, 'DELT2') 
    exp_MA_mean = data_RTSA.DELT2_mean;
    exp_MA_sd   = data_RTSA.DELT2_sd;
elseif strcmp(muscle_name, 'DELT3')
    exp_MA_mean = data_RTSA.DELT3_mean;
    exp_MA_sd   = data_RTSA.DELT3_sd;
end



%% Visualise after via point optimisation
figure(101);
scatter3(p_sim_0(1), p_sim_0(2), p_sim_0(3), 'o', 'filled','red')
scatter3(p_sim_opt(1), p_sim_opt(2), p_sim_opt(3), 'o', 'filled','black')

% Plot sphere to visualise constraint violation
theta = (0:0.01:1)*2*pi;
phi = (0:0.01:1)*pi;

[THETA,PHI]=meshgrid(theta,phi);
X1=radius.*cos(THETA).*sin(PHI) + p_sim_0(1);
Y1=radius.*sin(THETA).*sin(PHI) + p_sim_0(2);
Z1=radius.*cos(PHI) + p_sim_0(3);

figure(101);
surf(X1,Y1,Z1,...
    'FaceColor',[ 1 1 0],...
    'FaceAlpha', 0.2,...
    'EdgeColor', [0 0 0 ],...
    'EdgeAlpha', 0.1);


% Plot initial and optimised DELT1 MA

% Set via point to be the correct optimised position p_sim_opt
delt_via_downCast.set_location(Vec3(p_sim_opt(1), p_sim_opt(2), p_sim_opt(3)));
% Call ::Model.finalizeConnections() to....
osim_model.finalizeConnections();
% Re initialise the system to allow computing moment arm
new_state = osim_model.initSystem;

% % % % Re-calculate moment arms at each of the poses
% % % %%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 1 - 2.5 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % osim_model.updCoordinateSet().get('shoulder_elv').setValue(new_state, deg2rad(data_RTSA.angles(1)), true);
% % % osim_model.realizePosition(init_state);
% % % model_MA_optim.DELTx(1) = delt_GP.computeMomentArm(new_state,shoulder_elv);
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 2 - 30 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % osim_model.updCoordinateSet().get('shoulder_elv').setValue(new_state, deg2rad(data_RTSA.angles(2)), true);
% % % osim_model.realizePosition(new_state);
% % % model_MA_optim.DELTx(2) = delt_GP.computeMomentArm(new_state,shoulder_elv);
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 3 - 60 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % osim_model.updCoordinateSet().get('shoulder_elv').setValue(new_state, deg2rad(data_RTSA.angles(3)), true);
% % % osim_model.realizePosition(new_state);
% % % model_MA_optim.DELTx(3) = delt_GP.computeMomentArm(new_state,shoulder_elv);
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 4 - 90 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % osim_model.updCoordinateSet().get('shoulder_elv').setValue(new_state, deg2rad(data_RTSA.angles(4)), true);
% % % osim_model.realizePosition(new_state);
% % % model_MA_optim.DELTx(4) = delt_GP.computeMomentArm(new_state,shoulder_elv);
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 5 - 120 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % osim_model.updCoordinateSet().get('shoulder_elv').setValue(new_state, deg2rad(data_RTSA.angles(5)), true);
% % % osim_model.realizePosition(new_state);
% % % model_MA_optim.DELTx(5) = delt_GP.computeMomentArm(new_state,shoulder_elv);


figure(3);
hold on
title('deltoid moment arm for ith joint position (m)')
xlabel('Joint position');
ylabel('Moment arm');

hold on

for i_angle = 1:length(data_RTSA.angles)
   
    % Re-calculate moment arms at each of the poses
    %%%%%%%%%%%%%%%%%%% POS 1-5 - 2.5/30/60/90/120 DEG %%%%%%%%%%%%%%%%%%%%
    osim_model.updCoordinateSet().get('shoulder_elv').setValue(new_state, deg2rad(data_RTSA.angles(i_angle)), true);
    osim_model.realizePosition(init_state);
    model_MA_optim.DELTx(i_angle) = delt_GP.computeMomentArm(new_state,shoulder_elv);
    
    % Initial MA
    scatter(data_RTSA.angles(i_angle) , model_MA_init(i_angle),'filled','o','red');
    
    % Optimised MA    
    scatter(data_RTSA.angles(i_angle) , model_MA_optim.DELTx(i_angle),'filled','o','black');
    scatter(data_RTSA.angles(i_angle), exp_MA_mean(i_angle),'filled','o','magenta');
    scatter(data_RTSA.angles(i_angle) , exp_MA_mean(i_angle) + exp_MA_sd(i_angle),'+','magenta');
    scatter(data_RTSA.angles(i_angle) , exp_MA_mean(i_angle) - exp_MA_sd(i_angle), '+','magenta');

end

% % % scatter(data_RTSA.angles(1) , model_MA_init(1),'filled','o','red');
% % % scatter(data_RTSA.angles(2) , model_MA_init(2),'filled','o','red');
% % % scatter(data_RTSA.angles(3) , model_MA_init(3),'filled','o','red');
% % % scatter(data_RTSA.angles(4) , model_MA_init(4),'filled','o','red');
% % % scatter(data_RTSA.angles(5) , model_MA_init(5),'filled','o','red');
% % % 
% % % % Optimised MA
% % % %%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 1 - 2.5 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % scatter(data_RTSA.angles(1) , model_MA_optim.DELTx(1),'filled','o','black');
% % % 
% % % scatter(data_RTSA.angles(1), exp_MA_mean(1),'filled','o','magenta');
% % % scatter(data_RTSA.angles(1) , exp_MA_mean(1) + exp_MA_sd(1),'+','magenta');
% % % scatter(data_RTSA.angles(1) , exp_MA_mean(1) - exp_MA_sd(1), '+','magenta');
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 2 - 30 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % scatter(data_RTSA.angles(2) , model_MA_optim.DELTx(2),'filled','o','black');
% % % 
% % % scatter(data_RTSA.angles(2), exp_MA_mean(2),'filled','o','magenta');
% % % scatter(data_RTSA.angles(2) , exp_MA_mean(2) + exp_MA_sd(2),'+','magenta');
% % % scatter(data_RTSA.angles(2) , exp_MA_mean(2) - exp_MA_sd(2), '+','magenta');
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 3 - 60 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % scatter(data_RTSA.angles(3) , model_MA_optim.DELTx(3),'filled','o','black');
% % % 
% % % scatter(data_RTSA.angles(3), exp_MA_mean(3),'filled','o','magenta');
% % % scatter(data_RTSA.angles(3) , exp_MA_mean(3) + exp_MA_sd(3),'+','magenta');
% % % scatter(data_RTSA.angles(3) , exp_MA_mean(3) - exp_MA_sd(3), '+','magenta');
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 4 - 90 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % scatter(data_RTSA.angles(4) , model_MA_optim.DELTx(4),'filled','o','black');
% % % 
% % % scatter(data_RTSA.angles(4), exp_MA_mean(4),'filled','o','magenta');
% % % scatter(data_RTSA.angles(4) , exp_MA_mean(4) + exp_MA_sd(4),'+','magenta');
% % % scatter(data_RTSA.angles(4) , exp_MA_mean(4) - exp_MA_sd(4), '+','magenta');
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%% POSITION 5 - 120 DEG %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % scatter(data_RTSA.angles(5) , model_MA_optim.DELTx(5),'filled','o','black');
% % % 
% % % scatter(data_RTSA.angles(5), exp_MA_mean(5),'filled','o','magenta');
% % % scatter(data_RTSA.angles(5) , exp_MA_mean(5) + exp_MA_sd(5),'+','magenta');
% % % scatter(data_RTSA.angles(5) , exp_MA_mean(5) - exp_MA_sd(5), '+','magenta');

txt = ['\leftarrow J at optimum via point = ' num2str(J_opt)];
text(data_RTSA.angles(1)+0.01,model_MA_optim.DELTx(1),txt)

end