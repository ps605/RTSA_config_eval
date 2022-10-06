function plotViaOptResults(muscle_name, p_sim_0, p_sim_opt, delt_via_downCast, delt_GP, coord_for_MA, osim_model, data_RTSA, model_MA_init, J_opt, radius, flag_via_delete)
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

if flag_via_delete == false
    % Set via point to be the correct optimised position p_sim_opt
    delt_via_downCast.set_location(Vec3(p_sim_opt(1), p_sim_opt(2), p_sim_opt(3)));
    % Call ::Model.finalizeConnections() to....
    osim_model.finalizeConnections();
    % Re initialise the system to allow computing moment arm
    osim_model.initSystem;
end

figure;
hold on
title([muscle_name ' ' char(coord_for_MA.getName()) ' MA after via point opt'], 'Interpreter','none')
xlabel('Shoulder elevation angle (deg)');
ylabel('Moment arm (m)');
xlim([-5 125]);
xticks([2.5, 30, 60, 90, 120]);

hold on

for i_angle = 1:length(data_RTSA.angles)
   
    % Re-calculate moment arms at each of the poses
    %%%%%%%%%%%%%%%%%%% POS 1-5 - 2.5/30/60/90/120 DEG %%%%%%%%%%%%%%%%%%%%
    osim_model.updCoordinateSet().get('shoulder_elv').setValue(init_state, deg2rad(data_RTSA.angles(i_angle)), true);
    osim_model.realizePosition(init_state);
    model_MA_optim.DELTx(i_angle) = delt_GP.computeMomentArm(init_state,coord_for_MA);
    
    % Initial MA
    scatter(data_RTSA.angles(i_angle) , model_MA_init(i_angle),'filled','o','red');
    
    % Optimised MA    
    scatter(data_RTSA.angles(i_angle) , model_MA_optim.DELTx(i_angle),'filled','o','black');
    scatter(data_RTSA.angles(i_angle), exp_MA_mean(i_angle),'filled','o','magenta');
    scatter(data_RTSA.angles(i_angle) , exp_MA_mean(i_angle) + exp_MA_sd(i_angle),'+','magenta');
    scatter(data_RTSA.angles(i_angle) , exp_MA_mean(i_angle) - exp_MA_sd(i_angle), '+','magenta');

end


txt = ['\leftarrow J at optimum via point = ' num2str(J_opt) ' m for ' muscle_name];
text(data_RTSA.angles(1)+0.01,model_MA_optim.DELTx(1),txt)

end