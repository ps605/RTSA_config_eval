% Analysis of RTSA predictive simulations
% 1) JRF
% 2) Muscle Forces
% 3) Activations

clear
close all
clc
%% Setup


% Flags
flag_AnalysisTool   = false;
flag_NormTime       = true;
flag_HairTouch      = false; % O-RS_02
flag_LateralReach   = true; % O-RS_01 - sim_J46apu924ga
flag_UpwardReach    = false; % C-O-RS_01 - sim_J46apu924ga

% Import OpenSim 4.3 libraries
import org.opensim.modeling.*

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

jrf_colours = [255,109,106;
    240,230,140;
    222, 235, 247;
    221,160,221]./255;


jointF_to_plot = {'unrothum_on_scapula_in_scapula_offset_fx'...
    'unrothum_on_scapula_in_scapula_offset_fy',...
    'unrothum_on_scapula_in_scapula_offset_fz'};

muscleF_to_plot = {'DELT1',...
    'DELT2',...
    'DELT3'};

coords_to_plot = {
    'shoulder_elv'};

modes = {'m2',...
    'm4',...
    'm5',...
    'm6',...
    'm7',...
    'm9'};

SDs = {'-3',...
    '-1',...
    '1',...
    '3'};

if flag_HairTouch == true

    % Greens
    sd_colours = [229,245,224;
        161,217,155;
        65,171,93;
        0,109,44]./255;

    y_lims_force = {[-300, 200],...
        [-100, 500],...
        [0, 1200],...
        [0, 1200]};

    y_lims_muscle_force = {[0, 100],...
        [0, 600],...
        [0, 150]};

    y_lims_muscle_lengths = {[0.13, 0.22],...
        [0.10, 0.18],...
        [0.15, 0.19]};

    y_lims_muscle_moment = {%[0, 1], [-3, 0], [-2, 0];...
        [0, 5], [0, 10], [-1, 1];...
        %[-1.5, 0.2], [-7, 0], [-1.2, 0]...
        };

    y_lims_muscle_moment_arm = {[-5, 15], [-18, 0], [-30, 0];...
        [10, 50], [20, 50], [-15, 10];...
        [-20, 50], [-30, 0], [-15, 10]};

    task_name = 'HairTouch';

elseif flag_LateralReach == true

    % Reds
    sd_colours = [254,153,41;
        236,112,20;
        204,76,2;
        140,45,4]./255;

    y_lims_force = {[0, 500],...
        [0, 500],...
        [0, 500],...
        [0, 600]};

    y_lims_muscle_force = {[0, 100],...
        [0, 600],...
        [0, 150]};

    y_lims_muscle_lengths = {[0.13, 0.22],...
        [0.10, 0.18],...
        [0.15, 0.19]};

    y_lims_muscle_moment = {%[0, 1], [-3, 0], [-2, 0];...
        [0, 5], [0, 10], [-1, 1];...
        %[-1.5, 0.2], [-7, 0], [-1.2, 0]...
        };

    y_lims_muscle_moment_arm = {[10, 45], [10, 45], [-20, 5];...
        [10, 50], [20, 50], [-15, 10];...
        [-20, 50], [-30, 0], [-15, 10]};

    task_name = 'LateralReach';

elseif flag_UpwardReach == true

    % Blues
    sd_colours = [222, 235, 247;
        158, 202, 225;
        66, 146, 198;
        8, 69, 148]./255;

    y_lims_force = {[-400, 200],...
        [0, 500],...
        [0, 1000],...
        [0, 1100]};

    y_lims_muscle_force = {[0, 100],...
        [0, 600],...
        [0, 150]};

    y_lims_muscle_lengths = {[0.13, 0.22],...
        [0.10, 0.18],...
        [0.15, 0.19]};

    y_lims_muscle_moment = {%[0, 1], [-3, 0], [-2, 0];...
        [0, 5], [0, 10], [-1, 1];...
        %[-1.5, 0.2], [-7, 0], [-1.2, 0]...
        };

    y_lims_muscle_moment_arm = {[10, 35], [10, 35], [-20, 5];...
        [10, 50], [20, 50], [-15, 10];...
        [-20, 50], [-30, 0], [-15, 10]};

    task_name = 'UpwardReach';

end

analysis_folder = '../../OpenSim/Out/Moco/Analysis/ORS_01/';

% List simulation folders
sims = dir([analysis_folder '/sim_*']);

num_sims = numel(sims);

log_table = readtable('..\..\OpenSim\Out\Moco\Analysis\ORS_01\ORS_01_sims.txt');




%% Re-run AnalysisTool if needed
if flag_AnalysisTool == true

    for i_sim = 1:num_sims
        model_file = ['..\..\OpenSim\In\Models\RTSA_Adjusted\FSModel_GHJoint_' sims(i_sim).name(5:end) '.osim'];
        osim_model = Model();

        solution_file = [analysis_folder sims(i_sim).name '/MocoSol_' task_name '.sto'];

        states_storage=Storage(solution_file);

        % Set ArrayStr
        joints = ArrayStr();
        joints.set(0, 'unrothum');

        on_bodies = ArrayStr();
        on_bodies.set(0, 'parent');

        in_frames = ArrayStr();
        in_frames.set(0, '/jointset/unrothum/scapula_offset/');

        % Get time
        time_array = ArrayDouble;
        states_storage.getTimeColumn(time_array);

        % Set up AnalyzeTool
        analyzeTool=AnalyzeTool('..\..\OpenSim\In\Setup_files\Analysis\template_JRA_FR_MA.xml',0);
        analyzeTool.setName('Moco');
        analyzeTool.setInitialTime(0);
        analyzeTool.setFinalTime(time_array.getLast);
        %             analyzeTool.setStatesStorage(states_storage);
        analyzeTool.setStatesFileName(solution_file);
        %             analyzeTool.setModel(osim_model);
        analyzeTool.setModelFilename(model_file)
        analyzeTool.setResultsDir(['..\..\OpenSim\Out\Moco\Analysis\' sims(i_sim).name '\']);

        JR_Analysis = analyzeTool.updAnalysisSet.get(0);
        JR_downCast = JointReaction.safeDownCast(JR_Analysis);
        JR_downCast.setEndTime(time_array.getLast);
        JR_downCast.setInFrame(in_frames);
        JR_downCast.setOnBody(on_bodies);
        JR_downCast.setJointNames(joints);

        FR = analyzeTool.updAnalysisSet.get(1);
        FR.setEndTime(time_array.getLast);

        MA = analyzeTool.updAnalysisSet.get(2);
        MA.setEndTime(time_array.getLast);

        % Print and read back in as bug workaround
        analyzeTool.print('runAnalyzeTool.xml');

        runTool = AnalyzeTool('runAnalyzeTool.xml');
        runTool.run();

    end
end

%% Analysis Joint Reaction Analysis and kinematics
for i_sim = 1 : num_sims
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Import data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Import JRF data
    joint_reaction = importdata([analysis_folder sims(i_sim).name '/Moco_JointReaction_ReactionLoads.sto']);
    % Import States and Kinematics
    states = importdata([analysis_folder sims(i_sim).name '/MocoSol_' task_name '.sto']);
    shoulder_elv_theta.(sims(i_sim).name) = rad2deg(states.data(:,19));
    elv_angle_theta.(sims(i_sim).name) = rad2deg(states.data(:,18));
    shoulder_rot_theta.(sims(i_sim).name) = rad2deg(states.data(:,21));

    % Import Kinematics data
    %%% sim_JRF = importdata([analysis_folder sims(i_sim).name '/MocoSol_' task_name '.sto']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JRF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get time
    if flag_NormTime == true
        time_vec = linspace(0, 100, 101)';
        JRA.(sims(i_sim).name).time_og = 100*joint_reaction.data(:,1)/max(joint_reaction.data(:,1));
        JRA.(sims(i_sim).name).time = 100*joint_reaction.data(:,1)/max(joint_reaction.data(:,1));
    else
        JRA.(sims(i_sim).name).time = joint_reaction.data(:,1);
    end


    for i_joint = 1:numel(jointF_to_plot)
        % Get handle of coordinates to plot
        JRA.(sims(i_sim).name).pos(:, i_joint) = find(contains(joint_reaction.colheaders, jointF_to_plot{i_joint}));
        JRA.(sims(i_sim).name).label{1,i_joint} = joint_reaction.colheaders{JRA.(sims(i_sim).name).pos(i_joint)};

        if flag_NormTime == false
            if contains(JRA.(sims(i_sim).name).label{1,i_joint}, 'fz')
                JRA.(sims(i_sim).name).data(:, i_joint) = abs(joint_reaction.data(:,JRA.(sims(i_sim).name).pos(i_joint)));
            else
                JRA.(sims(i_sim).name).data(:, i_joint) = joint_reaction.data(:,JRA.(sims(i_sim).name).pos(i_joint));
            end
        elseif flag_NormTime == true

            if contains(JRA.(sims(i_sim).name).label{1,i_joint}, 'fz')
                JRA.(sims(i_sim).name).data(:, i_joint) = interp1(JRA.(sims(i_sim).name).time_og,...
                    abs(joint_reaction.data(:,JRA.(sims(i_sim).name).pos(i_joint))),...
                    time_vec);
            else
                JRA.(sims(i_sim).name).data(:, i_joint) = interp1(JRA.(sims(i_sim).name).time_og,...
                    joint_reaction.data(:,JRA.(sims(i_sim).name).pos(i_joint)),...
                    time_vec);
            end

            JRA.(sims(i_sim).name).time = time_vec;

            
        end

        % Create JRF plots
        figure (i_joint)

        plot(JRA.(sims(i_sim).name).time, JRA.(sims(i_sim).name).data(:, i_joint), 'LineWidth', 1.5, 'Color', jrf_colours(i_joint,:))

        hold on;
        xlabel('Movement Duration (%)');
        ylabel(['JRF [' jointF_to_plot{i_joint}(end) '] (N)']);
        ylim(y_lims_force{i_joint});
        hold on;



    end

    % Normalise Kinematics
    shoulder_elv_theta.(sims(i_sim).name) = interp1(JRA.(sims(i_sim).name).time_og,...
        shoulder_elv_theta.(sims(i_sim).name),...
        time_vec);

    elv_angle_theta.(sims(i_sim).name) = interp1(JRA.(sims(i_sim).name).time_og,...
        elv_angle_theta.(sims(i_sim).name),...
        time_vec);

    shoulder_rot_theta.(sims(i_sim).name) = interp1(JRA.(sims(i_sim).name).time_og,...
        shoulder_rot_theta.(sims(i_sim).name),...
        time_vec);

    % Calculate resultant Force
    JRA.(sims(i_sim).name).F_res = sqrt(JRA.(sims(i_sim).name).data(:, 1).^2 + JRA.(sims(i_sim).name).data(:, 2).^2 + JRA.(sims(i_sim).name).data(:, 3).^2);

    figure(4)

    plot(JRA.(sims(i_sim).name).time, JRA.(sims(i_sim).name).F_res, 'LineWidth', 1.5, 'Color', jrf_colours(4,:))

    hold on;
    %scatter(JRA.(sims(i_sim).name).time(max_point + round(numel(JRA.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
    xlabel('Movement Duration (%)');
    ylabel('JRF [Resultant] (N)');
    ylim(y_lims_force{4});
    hold on;

    % Plot Joint Kinematics

    figure (11)

    plot(JRA.(sims(i_sim).name).time, shoulder_elv_theta.(sims(i_sim).name), 'LineWidth', 1.5, 'Color', jrf_colours(1,:))

    hold on;
    xlabel('Movement Duration (%)');
    ylabel(['Glenohumeral Elevation (deg)']);
    hold on;

    figure (12)

    plot(JRA.(sims(i_sim).name).time, elv_angle_theta.(sims(i_sim).name), 'LineWidth', 1.5, 'Color', jrf_colours(2,:))

    hold on;
    xlabel('Movement Duration (%)');
    ylabel(['Elevation angle (deg)']);
    hold on;

    figure (13)

    plot(JRA.(sims(i_sim).name).time, shoulder_rot_theta.(sims(i_sim).name), 'LineWidth', 1.5, 'Color', jrf_colours(3,:))

    hold on;
    xlabel('Movement Duration (%)');
    ylabel(['Shoulder Rotation (deg)']);
    hold on;

end

% Plot mean scapula data on top of other trials
figure(1)
plot(JRA.sim_J46apu924ga.time, JRA.sim_J46apu924ga.data(:, 1), 'LineWidth', 1.5, 'Color', 'black')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Fx_norm.tiff','tiff')

figure(2)
plot(JRA.sim_J46apu924ga.time, JRA.sim_J46apu924ga.data(:, 2), 'LineWidth', 1.5, 'Color', 'black')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Fy_norm.tiff','tiff')

figure(3)
plot(JRA.sim_J46apu924ga.time, JRA.sim_J46apu924ga.data(:, 3), 'LineWidth', 1.5, 'Color', 'black')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Fz_norm.tiff','tiff')

figure(4)
plot(JRA.sim_J46apu924ga.time, JRA.sim_J46apu924ga.F_res(:, 1), 'LineWidth', 1.5, 'Color', 'black')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(4), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Fres_norm.tiff','tiff')

figure(11)
plot(JRA.sim_J46apu924ga.time, shoulder_elv_theta.sim_J46apu924ga, 'LineWidth', 1.5, 'Color', 'black')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(11), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\shoulder_elv_norm.tiff','tiff')

figure(12)
plot(JRA.sim_J46apu924ga.time, elv_angle_theta.sim_J46apu924ga, 'LineWidth', 1.5, 'Color', 'black')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(12), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\elv_ang_norm.tiff','tiff')

figure(13)
plot(JRA.sim_J46apu924ga.time, shoulder_rot_theta.sim_J46apu924ga, 'LineWidth', 1.5, 'Color', 'black')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(13), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\shoulder_rot_norm.tiff','tiff')

close all

%% By mode
for i_mode = 1:numel(modes)

    for i_sd = 1:numel(SDs)

        morphology = [modes{i_mode} '_' SDs{i_sd}];

        idx_sim_table =  find(contains(log_table.Scapula_morphology, morphology));

        for i_joint = 1:3

            figure (i_joint)
            subplot(2,3, i_mode)
            plot(JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).time,...
                JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_joint), ...
                'LineWidth', 1.5,...
                'Color', sd_colours(i_sd,:))

            %             plot(shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_sim_table}]),...
            %                 JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_joint), ...
            %                 'LineWidth', 1.5,...
            %                 'Color', sd_colours(i_sd,:))

            hold on;
            %             scatter(JRA.(sims(i_sim).name).time(max_point + round(numel(JRA.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
            xlabel('Movement Duration (%)');
            ylabel(['JRF [' jointF_to_plot{i_joint}(end) '] (N)']);
            ylim(y_lims_force{i_joint});
            title(modes{i_mode})
            hold on;

        end

        figure(4)
        subplot(2,3, i_mode)
        plot(JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).time, ...
            JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).F_res ,...
            'LineWidth', 1.5, ...
            'Color', sd_colours(i_sd,:))


        %         plot(shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_sim_table}]), ...
        %             JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).F_res ,...
        %             'LineWidth', 1.5, ...
        %             'Color', sd_colours(i_sd,:))

        hold on;
        %scatter(JRA.(sims(i_sim).name).time(max_point + round(numel(JRA.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
        xlabel('Movement Duration (%)','FontWeight', 'bold');
        ylabel('JRF [Resultant] (N)', 'FontWeight', 'bold');
        ylim(y_lims_force{4});
        %xlim([0 80])
        title(modes{i_mode})
        hold on;

        % Plot Force ratio
        figure(5)

        fr = sqrt(JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, 3).^2)./...
            (sqrt(JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, 1).^2) + sqrt(JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, 2).^2));


        subplot(2,3, i_mode)
        plot(JRA.(sims(i_sim).name).time, fr, 'LineWidth', 1.5, 'Color', sd_colours(i_sd,:),'LineStyle', '--')
        hold on
        ylabel('JRF Compressive to Shear Ratio', 'FontWeight', 'bold');
        xlabel('Movement Duration (%)','FontWeight', 'bold');
        ylim([0, 2])
        title(modes{i_mode})


    end

    % Average morphology data
    figure(1)
    plot(JRA.sim_J46apu924ga.time, JRA.sim_J46apu924ga.data(:, 1), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
    figure(2)
    plot(JRA.sim_J46apu924ga.time, JRA.sim_J46apu924ga.data(:, 2), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
    figure(3)
    plot(JRA.sim_J46apu924ga.time, JRA.sim_J46apu924ga.data(:, 3), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
    figure(4)
    plot(JRA.sim_J46apu924ga.time, JRA.sim_J46apu924ga.F_res(:, 1), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')

    figure(5)
    fr_av = sqrt(JRA.sim_J46apu924ga.data(:, 3).^2)./...
        (sqrt(JRA.sim_J46apu924ga.data(:, 1).^2) + sqrt(JRA.sim_J46apu924ga.data(:, 2).^2));
    plot(JRA.sim_J46apu924ga.time, fr_av, 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
end

figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Fx_vs_NORM_modes.tiff','tiff')
figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Fy_vs_NORM_modes.tiff','tiff')
figure(3)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Fz_vs_NORM_modes.tiff','tiff')
figure(4)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(4), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Fres_vs_NORM_modes.tiff','tiff')
figure(5)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(5), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\FR_vs_NORM_modes.tiff','tiff')



close all
%% Analysis muscle forces
for i_sim = 1: num_sims

    % Import ForceReporter data
    tendon_force = importdata([analysis_folder sims(i_sim).name '/Moco_MuscleAnalysis_TendonForce.sto']);

    % Get time
    TF.(sims(i_sim).name).time = tendon_force.data(:,1);

    for i_mus = 1:numel(muscleF_to_plot)
        % Get handle of coordinates to plot
        TF.(sims(i_sim).name).pos(:, i_mus) = find(contains(tendon_force.colheaders, muscleF_to_plot{i_mus}));
        TF.(sims(i_sim).name).label{1,i_mus} = tendon_force.colheaders{TF.(sims(i_sim).name).pos(i_mus)};

        TF.(sims(i_sim).name).data(:, i_mus) = tendon_force.data(:,TF.(sims(i_sim).name).pos(i_mus));


        % Create q plots
        figure (i_mus)


        plot(TF.(sims(i_sim).name).time, TF.(sims(i_sim).name).data(:, i_mus), 'LineWidth', 1.5, 'Color', jrf_colours(i_mus,:))

        hold on;
        %         scatter(TF.(sims(i_sim).name).time(max_point + round(numel(TF.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
        xlabel('Time (s)');
        ylabel(['Total Muscle (Tendon) Force [' muscleF_to_plot{i_mus} '] (N)']);
        hold on;
    end

end


% Plot mean scapula data on top of other trials
figure(1)
plot(TF.sim_J46apu924ga.time, TF.sim_J46apu924ga.data(:, 1), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Force_DELT1.tiff','tiff')
figure(2)
plot(TF.sim_J46apu924ga.time, TF.sim_J46apu924ga.data(:, 2), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Force_DELT2.tiff','tiff')
figure(3)
plot(TF.sim_J46apu924ga.time, TF.sim_J46apu924ga.data(:, 3), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Force_DELT3.tiff','tiff')

close all

%% By mode
for i_mode = 1:numel(modes)

    for i_sd = 1:numel(SDs)

        morphology = [modes{i_mode} '_' SDs{i_sd}];

        idx_sim_table =  find(contains(log_table.Scapula_morphology, morphology));

        for i_mus = 1:3

            figure (i_mus)
            subplot(2,3, i_mode)
            %             plot(TF.(['sim_' log_table.Model_Hash{idx_sim_table}]).time,...
            %                 TF.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus), ...
            %                 'LineWidth', 1.5,...
            %                 'Color', sd_colours(i_sd,:))

            plot(shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_sim_table}]),...
                TF.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus), ...
                'LineWidth', 1.5,...
                'Color', sd_colours(i_sd,:))

            hold on;
            %             scatter(JRA.(sims(i_sim).name).time(max_point + round(numel(JRA.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
            xlabel('Glenohumeral ABD (deg)');
            ylabel(['Total Muscle (Tendon) Force [' muscleF_to_plot{i_mus} '] (N)']);
            ylim(y_lims_muscle_force{i_mus});
            title(modes{i_mode})
            hold on;

        end


    end

    figure(1)
    plot(shoulder_elv_theta.sim_J46apu924ga, TF.sim_J46apu924ga.data(:, 1), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
    figure(2)
    plot(shoulder_elv_theta.sim_J46apu924ga, TF.sim_J46apu924ga.data(:, 2), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
    figure(3)
    plot(shoulder_elv_theta.sim_J46apu924ga, TF.sim_J46apu924ga.data(:, 3), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')

end


figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Force_vs_ABD_DELT1_modes.tiff','tiff')
figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Force_vs_ABD_DELT2_modes.tiff','tiff')
figure(3)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Force_vs_ABD_DELT3_modes.tiff','tiff')

close all

%% Muscle Lengths
for i_sim = 1: num_sims

    % Import ForceReporter data
    muscle_length = importdata([analysis_folder sims(i_sim).name '/Moco_MuscleAnalysis_Length.sto']);

    % Get time
    ML.(sims(i_sim).name).time = muscle_length.data(:,1);

    for i_mus = 1:numel(muscleF_to_plot)
        % Get handle of coordinates to plot
        ML.(sims(i_sim).name).pos(:, i_mus) = find(contains(muscle_length.colheaders, muscleF_to_plot{i_mus}));
        ML.(sims(i_sim).name).label{1,i_mus} = muscle_length.colheaders{ML.(sims(i_sim).name).pos(i_mus)};

        ML.(sims(i_sim).name).data(:, i_mus) = muscle_length.data(:,ML.(sims(i_sim).name).pos(i_mus));


        % Create q plots
        figure (i_mus)


        plot(shoulder_elv_theta.(sims(i_sim).name), ML.(sims(i_sim).name).data(:, i_mus), 'LineWidth', 1.5, 'Color', jrf_colours(i_mus,:))

        hold on;
        %         scatter(TF.(sims(i_sim).name).time(max_point + round(numel(TF.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
        xlabel('Time (s)');
        ylabel([' Muscle Length [' muscleF_to_plot{i_mus} '] (m)']);
        hold on;
    end

end


% Plot mean scapula data on top of other trials
figure(1)
plot(ML.sim_J46apu924ga.time, ML.sim_J46apu924ga.data(:, 1), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Length_DELT1.tiff','tiff')
figure(2)
plot(ML.sim_J46apu924ga.time, ML.sim_J46apu924ga.data(:, 2), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Length_DELT2.tiff','tiff')
figure(3)
plot(ML.sim_J46apu924ga.time, ML.sim_J46apu924ga.data(:, 3), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Length_DELT3.tiff','tiff')

close all

%% By mode

for i_mode = 1:numel(modes)

    for i_sd = 1:numel(SDs)

        morphology = [modes{i_mode} '_' SDs{i_sd}];

        idx_sim_table =  find(contains(log_table.Scapula_morphology, morphology));

        for i_mus = 1:3

            figure (i_mus)
            subplot(2,3, i_mode)
            %             plot(ML.(['sim_' log_table.Model_Hash{idx_sim_table}]).time,...
            %                 ML.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus), ...
            %                 'LineWidth', 1.5,...
            %                 'Color', sd_colours(i_sd,:))

            plot(shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_sim_table}]),...
                ML.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus), ...
                'LineWidth', 1.5,...
                'Color', sd_colours(i_sd,:))

            hold on;
            %             scatter(JRA.(sims(i_sim).name).time(max_point + round(numel(JRA.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
            xlabel('Glenohumeral ABD (deg)');
            ylabel(['Muscle Length [' muscleF_to_plot{i_mus} '] (m)']);
            ylim(y_lims_muscle_lengths{i_mus});
            title(modes{i_mode})
            hold on;

        end


    end

    figure(1)
    plot(shoulder_elv_theta.sim_J46apu924ga, ML.sim_J46apu924ga.data(:, 1), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
    figure(2)
    plot(shoulder_elv_theta.sim_J46apu924ga, ML.sim_J46apu924ga.data(:, 2), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
    figure(3)
    plot(shoulder_elv_theta.sim_J46apu924ga, ML.sim_J46apu924ga.data(:, 3), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')

end


figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Length_vs_ABD_DELT1_modes.tiff','tiff')
figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Length_vs_ABD_DELT2_modes.tiff','tiff')
figure(3)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Length_vs_ABD_DELT3_modes.tiff','tiff')
close all

%% Muscle Activations
for i_sim = 1: num_sims

    % Import ForceReporter data
    muscle_activation = importdata([analysis_folder sims(i_sim).name '/MocoSol_' task_name '.sto']);

    % Get time
    MAct.(sims(i_sim).name).time = muscle_activation.data(:,1);

    for i_mus = 1:numel(muscleF_to_plot)
        % Get handle of coordinates to plot
        MAct.(sims(i_sim).name).pos(:, i_mus) = find(contains(muscle_activation.colheaders, [muscleF_to_plot{i_mus} '/activation']));
        MAct.(sims(i_sim).name).label{1,i_mus} = muscle_activation.colheaders{MAct.(sims(i_sim).name).pos(i_mus)};

        MAct.(sims(i_sim).name).data(:, i_mus) = muscle_activation.data(:,MAct.(sims(i_sim).name).pos(i_mus));


        % Create q plots
        figure (i_mus)


        plot(MAct.(sims(i_sim).name).time, MAct.(sims(i_sim).name).data(:, i_mus), 'LineWidth', 1.5, 'Color', jrf_colours(i_mus,:))

        hold on;
        %         scatter(TF.(sims(i_sim).name).time(max_point + round(numel(TF.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
        xlabel('Time (s)');
        ylabel([' Muscle Activation [' muscleF_to_plot{i_mus} ']']);
        ylim([0, 1]);
        hold on;
    end

end


% Plot mean scapula data on top of other trials
figure(1)
plot(MAct.sim_J46apu924ga.time, MAct.sim_J46apu924ga.data(:, 1), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Activation_DELT1.tiff','tiff')
figure(2)
plot(MAct.sim_J46apu924ga.time, MAct.sim_J46apu924ga.data(:, 2), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Activation_DELT2.tiff','tiff')
figure(3)
plot(MAct.sim_J46apu924ga.time, MAct.sim_J46apu924ga.data(:, 3), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Activation_DELT3.tiff','tiff')

close all

%% By mode

for i_mode = 1:numel(modes)

    for i_sd = 1:numel(SDs)

        morphology = [modes{i_mode} '_' SDs{i_sd}];

        idx_sim_table =  find(contains(log_table.Scapula_morphology, morphology));

        for i_mus = 1:3

            figure (i_mus)
            subplot(2,3, i_mode)
            plot(MAct.(['sim_' log_table.Model_Hash{idx_sim_table}]).time,...
                MAct.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus), ...
                'LineWidth', 1.5,...
                'Color', sd_colours(i_sd,:))

            hold on;
            %             scatter(JRA.(sims(i_sim).name).time(max_point + round(numel(JRA.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
            xlabel('Time (s)');
            ylabel(['Muscle Activation [' muscleF_to_plot{i_mus} ']']);
            ylim([0, 1]);
            title(modes{i_mode})
            hold on;

        end


    end

    figure(1)
    plot(MAct.sim_J46apu924ga.time, MAct.sim_J46apu924ga.data(:, 1), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
    figure(2)
    plot(MAct.sim_J46apu924ga.time, MAct.sim_J46apu924ga.data(:, 2), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
    figure(3)
    plot(MAct.sim_J46apu924ga.time, MAct.sim_J46apu924ga.data(:, 3), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')

end


figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Activation_DELT1_modes.tiff','tiff')
figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Activation_DELT2_modes.tiff','tiff')
figure(3)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), '..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Activation_DELT3_modes.tiff','tiff')

%% Analysis muscle moment
for i_coord = 1:numel(coords_to_plot)

    coord = coords_to_plot{i_coord};

    for i_sim = 1: num_sims

        % Import ForceReporter data
        muscle_moment = importdata([analysis_folder sims(i_sim).name '/Moco_MuscleAnalysis_Moment_' coord '.sto']);

        % Get time
        MM.(sims(i_sim).name).time = muscle_moment.data(:,1);

        for i_mus = 1:numel(muscleF_to_plot)
            % Get handle of coordinates to plot
            MM.(sims(i_sim).name).pos(:, i_mus) = find(contains(muscle_moment.colheaders, muscleF_to_plot{i_mus}));
            MM.(sims(i_sim).name).label{1,i_mus} = muscle_moment.colheaders{MM.(sims(i_sim).name).pos(i_mus)};

            MM.(sims(i_sim).name).data(:, i_mus) = muscle_moment.data(:,MM.(sims(i_sim).name).pos(i_mus));


            % Create q plots
            figure (i_mus)


            plot(shoulder_elv_theta.(sims(i_sim).name), MM.(sims(i_sim).name).data(:, i_mus), 'LineWidth', 1.5, 'Color', jrf_colours(i_mus,:))

            hold on;
            %         scatter(TF.(sims(i_sim).name).time(max_point + round(numel(TF.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
            xlabel('Time (s)');
            ylabel(['Muscle Moment [' muscleF_to_plot{i_mus} '] (Nm)']);
            hold on;
        end

    end


    % Plot mean scapula data on top of other trials
    figure(1)
    plot(MM.sim_J46apu924ga.time, MM.sim_J46apu924ga.data(:, 1), 'LineWidth', 1.5, 'Color', 'black')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(1), ['..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Moment_DELT1_' coord '.tiff'],'tiff')
    figure(2)
    plot(MM.sim_J46apu924ga.time, MM.sim_J46apu924ga.data(:, 2), 'LineWidth', 1.5, 'Color', 'black')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(2), ['..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Moment_DELT2_' coord '.tiff'],'tiff')
    figure(3)
    plot(MM.sim_J46apu924ga.time, MM.sim_J46apu924ga.data(:, 3), 'LineWidth', 1.5, 'Color', 'black')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(3), ['..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Moment_DELT3_' coord '.tiff'],'tiff')

    close all

    %% By mode
    for i_mode = 1:numel(modes)

        for i_sd = 1:numel(SDs)

            morphology = [modes{i_mode} '_' SDs{i_sd}];

            idx_sim_table =  find(contains(log_table.Scapula_morphology, morphology));

            for i_mus = 1:3

                figure (i_mus)
                subplot(2,3, i_mode)
                %                 plot(MM.(['sim_' log_table.Model_Hash{idx_sim_table}]).time,...
                %                     MM.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus), ...
                %                     'LineWidth', 1.5,...
                %                     'Color', sd_colours(i_sd,:))

                plot(shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_sim_table}]),...
                    MM.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus), ...
                    'LineWidth', 1.5,...
                    'Color', sd_colours(i_sd,:))

                hold on;
                %             scatter(JRA.(sims(i_sim).name).time(max_point + round(numel(JRA.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
                xlabel('Glenohumeral ABD (deg)');
                ylabel(['Muscle Moment [' muscleF_to_plot{i_mus} '] (Nm)']);
                ylim(y_lims_muscle_moment{i_coord, i_mus});
                title(modes{i_mode})
                hold on;

            end


        end

        figure(1)
        plot(shoulder_elv_theta.sim_J46apu924ga, MM.sim_J46apu924ga.data(:, 1), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
        figure(2)
        plot(shoulder_elv_theta.sim_J46apu924ga, MM.sim_J46apu924ga.data(:, 2), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
        figure(3)
        plot(shoulder_elv_theta.sim_J46apu924ga, MM.sim_J46apu924ga.data(:, 3), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')

    end


    figure(1)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(1), ['..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Moment_vs_ABD_DELT1_' coord '_modes.tiff'],'tiff')
    figure(2)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(2), ['..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Moment_vs_ABD_DELT2_' coord '_modes.tiff'],'tiff')
    figure(3)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(3), ['..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\Moment_vs_ABD_DELT3_' coord '_modes.tiff'],'tiff')

    close all
end

%% Analysis muscle moment arms
for i_coord = 1:numel(coords_to_plot)

    coord = coords_to_plot{i_coord};

    for i_sim = 1: num_sims

        % Import ForceReporter data
        muscle_moment_arm = importdata([analysis_folder sims(i_sim).name '/Moco_MuscleAnalysis_MomentArm_' coord '.sto']);

        % Get time
        MA.(sims(i_sim).name).time = muscle_moment_arm.data(:,1);

        for i_mus = 1:numel(muscleF_to_plot)
            % Get handle of coordinates to plot
            MA.(sims(i_sim).name).pos(:, i_mus) = find(contains(muscle_moment_arm.colheaders, muscleF_to_plot{i_mus}));
            MA.(sims(i_sim).name).label{1,i_mus} = muscle_moment_arm.colheaders{MA.(sims(i_sim).name).pos(i_mus)};

            MA.(sims(i_sim).name).data(:, i_mus) = muscle_moment_arm.data(:,MA.(sims(i_sim).name).pos(i_mus));


            % Create q plots
            figure (i_mus)


            plot(shoulder_elv_theta.(sims(i_sim).name),...
                MA.(sims(i_sim).name).data(:, i_mus)*1000, ...
                'LineWidth', 1.5, ...
                'Color', jrf_colours(i_mus,:))

            hold on;
            %         scatter(TF.(sims(i_sim).name).time(max_point + round(numel(TF.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
            xlabel('Glenohumeral ABD (deg)');
            ylabel(['Muscle Moment Arm [' muscleF_to_plot{i_mus} '] (mm)']);
            hold on;
        end

    end


    % Plot mean scapula data on top of other trials
    figure(1)
    plot(shoulder_elv_theta.sim_J46apu924ga, MA.sim_J46apu924ga.data(:, 1)*1000, 'LineWidth', 1.5, 'Color', 'black')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(1), ['..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\MomentArm_DELT1_' coord '.tiff'],'tiff')
    figure(2)
    plot(shoulder_elv_theta.sim_J46apu924ga, MA.sim_J46apu924ga.data(:, 2)*1000, 'LineWidth', 1.5, 'Color', 'black')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(2), ['..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\MomentArm_DELT2_' coord '.tiff'],'tiff')
    figure(3)
    plot(shoulder_elv_theta.sim_J46apu924ga, MA.sim_J46apu924ga.data(:, 3)*1000, 'LineWidth', 1.5, 'Color', 'black')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(3), ['..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\MomentArm_DELT3_' coord '.tiff'],'tiff')

    close all

    %% By mode
    for i_mode = 1:numel(modes)

        for i_sd = 1:numel(SDs)

            morphology = [modes{i_mode} '_' SDs{i_sd}];

            idx_sim_table =  find(contains(log_table.Scapula_morphology, morphology));

            for i_mus = 1:3

                figure (i_mus)
                subplot(2,3, i_mode)
                plot(shoulder_elv_theta.(sims(i_sim).name),...
                    MA.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus)*1000, ...
                    'LineWidth', 1.5,...
                    'Color', sd_colours(i_sd,:))

                hold on;
                %             scatter(JRA.(sims(i_sim).name).time(max_point + round(numel(JRA.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
                xlabel('Glenohumeral ABD (deg)');
                ylabel(['Muscle Moment Arm [' muscleF_to_plot{i_mus} '] (mm)']);
                ylim(y_lims_muscle_moment_arm{i_coord, i_mus});
                title(modes{i_mode})
                hold on;

            end


        end

        figure(1)
        plot(shoulder_elv_theta.sim_J46apu924ga, MA.sim_J46apu924ga.data(:, 1)*1000, 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
        figure(2)
        plot(shoulder_elv_theta.sim_J46apu924ga, MA.sim_J46apu924ga.data(:, 2)*1000, 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
        figure(3)
        plot(shoulder_elv_theta.sim_J46apu924ga, MA.sim_J46apu924ga.data(:, 3)*1000, 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')

    end


    figure(1)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(1), ['..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\MomentArm_DELT1_' coord '_modes.tiff'],'tiff')
    figure(2)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(2), ['..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\MomentArm_DELT2_' coord '_modes.tiff'],'tiff')
    figure(3)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(3), ['..\..\OpenSim\Out\Moco\Analysis\Plots\ORS_01\MomentArm_DELT3_' coord '_modes.tiff'],'tiff')

    close all
end