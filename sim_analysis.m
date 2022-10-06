% Analysis of RTSA predictive simulations
% 1) JRF
% 2) Muscle Forces
% 3) Activations

clear
close all
clc
%% Setup
muscles_to_plot = {'SUPSP',...
    'INFSP',...
    'SUBSC',...
    'TMIN',...
    'TRP1',...
    'TRP2',...
    'TRP3',...
    'DELT1',...
    'DELT2',...
    'DELT3',...
    'LVS',...
    'PECM1',...
    'PECM2',...
    'PECM3'};

muscle_colours = [215 181 216;
    223 101 176;
    221 28 119;
    152 0 67;
    178 226 226;
    102 194 164;
    44 162 95;
    254 204 92;
    253 141 60;
    240 59 32;
    117 107 177;
    65 182 196;
    44 127 184;
    37 52 148]./255;

jrf_colours = [240 59 32;
    254 217 118
    44 127 184;
    221 28 119]./255;

jointF_to_plot = {'shoulder0_on_scapphant_in_glenoid_centre_fx'...
    'shoulder0_on_scapphant_in_glenoid_centre_fy',...
    'shoulder0_on_scapphant_in_glenoid_centre_fz'};

analysis_folder = '../../OpenSim/Out/Moco/Analysis/';

% List simulation folders
sims = dir([analysis_folder '/sim_*']);

num_sims = numel(sims);

%% Analysis
for i_sim = 1 : num_sims
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Import data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Import JRF data
    joint_reaction = importdata([analysis_folder sims(i_sim).name '/Moco_JointReaction_ReactionLoads.sto']);

    % Import ForceReporter data
    force_reporter = importdata([analysis_folder sims(i_sim).name '/Moco_ForceReporter_forces.sto']);

    % Import Kinematics data
    %%% sim_JRF = importdata([analysis_folder sims(i_sim).name '/MocoSol_LateralReach.sto']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JRF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get time
    JRA.time = joint_reaction.data(:,1);
    
    for i_joint = 1:numel(jointF_to_plot)
        % Get handle of coordinates to plot
        JRA.pos(:, i_joint) = find(contains(joint_reaction.colheaders, jointF_to_plot{i_joint}));
        JRA.label{1,i_joint} = joint_reaction.colheaders{JRA.pos(i_joint)};

        if contains(JRA.label{1,i_joint}, 'fz')
            JRA.data(:, i_joint) = abs(joint_reaction.data(:,JRA.pos(i_joint)));
        else
            JRA.data(:, i_joint) = joint_reaction.data(:,JRA.pos(i_joint));
        end
        
        % Get max force values
        [max_F_final(i_sim, i_joint), max_point] = max(JRA.data(round(numel(JRA.data(:, i_joint))*0.6):end,i_joint));

        % Create q plots
        figure (i_joint)
        plot(JRA.time, JRA.data(:, i_joint), 'LineWidth', 1.5, 'Color', jrf_colours(i_joint,:))
        hold on;
        scatter(JRA.time(max_point + round(numel(JRA.data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
        xlabel('Time (s)');
        ylabel(['JRF [' jointF_to_plot{i_joint}(end) '] (N)']);
        hold on;

    end


end
