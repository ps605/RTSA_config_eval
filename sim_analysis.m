% Analysis of RTSA predictive simulations
% 1) JRF
% 2) Muscle Forces
% 3) Activations

clear
close all
clc
%% Setup

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

jrf_colours = [240 59 32;
    254 217 118
    44 127 184;
    221 28 119]./255;
% Blues
sd_colours = [222, 235, 247;
    158, 202, 225;
    66, 146, 198;
    8, 69, 148]./255;
% Reds
% sd_colours = [254,153,41;
%     236,112,20;
%     204,76,2;
%     140,45,4]./255;

jointF_to_plot = {'unrothum_on_scapula_in_scapula_offset_fx'...
    'unrothum_on_scapula_in_scapula_offset_fy',...
    'unrothum_on_scapula_in_scapula_offset_fz'};

muscleF_to_plot = {'DELT1',...
    'DELT2',...
    'DELT3'};

coords_to_plot = {
    'shoulder_elv'};

modes = {%'m2',...
    %'m4',...
    %'m5',...
    'm6',...
    'm7',...
    'm9'};

SDs = {'-3',...
    '-1',...
    '1',...
    '3'};

y_lims_force = {[-300, 200],...
    [0, 500],...
    [0, 600],...
    [0, 600]};

y_lims_muscle_force = {[0, 100],...
    [0, 400],...
    [0, 200]};

y_lims_muscle_lengths = {[0.15, 0.22],...
    [0.11, 0.18],...
    [0.15, 0.18]};

y_lims_muscle_moment = {%[0, 1], [-3, 0], [-2, 0];...
    [0, 5], [0, 10], [-1, 1];...
    %[-1.5, 0.2], [-7, 0], [-1.2, 0]...
    };

y_lims_muscle_moment_arm = {[-5, 15], [-18, 0], [-30, 0];...
    [10, 50], [20, 50], [-15, 10];...
    [-20, 50], [-30, 0], [-15, 10]};

analysis_folder = '../../OpenSim/Out/Moco/Analysis/CORS_01/';

% List simulation folders
sims = dir([analysis_folder '/sim_*']);

num_sims = numel(sims);

log_table = readtable('..\..\OpenSim\Out\Moco\Analysis\CORS_01\CORS_01_sims.xlsx');

task_name = 'LateralReach';

% Flags
flag_AnalysisTool = false;

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

%% Analysis Joint Reaction Analysis
for i_sim = 1 : num_sims
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Import data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Import JRF data
    joint_reaction = importdata([analysis_folder sims(i_sim).name '/Moco_JointReaction_ReactionLoads.sto']);
    states = importdata([analysis_folder sims(i_sim).name '/MocoSol_UpwardReach.sto']);
    shoulder_elv_theta.(sims(i_sim).name) = rad2deg(states.data(:,19));

    % Import Kinematics data
    %%% sim_JRF = importdata([analysis_folder sims(i_sim).name '/MocoSol_LateralReach.sto']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JRF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get time
    JRA.(sims(i_sim).name).time = joint_reaction.data(:,1);
    
    for i_joint = 1:numel(jointF_to_plot)
        % Get handle of coordinates to plot
        JRA.(sims(i_sim).name).pos(:, i_joint) = find(contains(joint_reaction.colheaders, jointF_to_plot{i_joint}));
        JRA.(sims(i_sim).name).label{1,i_joint} = joint_reaction.colheaders{JRA.(sims(i_sim).name).pos(i_joint)};

        if contains(JRA.(sims(i_sim).name).label{1,i_joint}, 'fz')
            JRA.(sims(i_sim).name).data(:, i_joint) = abs(joint_reaction.data(:,JRA.(sims(i_sim).name).pos(i_joint)));
        else
            JRA.(sims(i_sim).name).data(:, i_joint) = joint_reaction.data(:,JRA.(sims(i_sim).name).pos(i_joint));
        end

        % Get max force values
        [max_F_final(i_sim, i_joint), max_point] = max(JRA.(sims(i_sim).name).data(round(numel(JRA.(sims(i_sim).name).data(:, i_joint))*0.6):end,i_joint));

        % Create q plots
        figure (i_joint)


        plot(JRA.(sims(i_sim).name).time, JRA.(sims(i_sim).name).data(:, i_joint), 'LineWidth', 1.5, 'Color', jrf_colours(i_joint,:))

        hold on;
        scatter(JRA.(sims(i_sim).name).time(max_point + round(numel(JRA.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
        xlabel('Time (s)');
        ylabel(['JRF [' jointF_to_plot{i_joint}(end) '] (N)']);
        ylim(y_lims_force{i_joint});
        hold on;

    end

    JRA.(sims(i_sim).name).F_res = sqrt(JRA.(sims(i_sim).name).data(:, 1).^2 + JRA.(sims(i_sim).name).data(:, 2).^2 + JRA.(sims(i_sim).name).data(:, 3).^2);
    figure(4)

    plot(JRA.(sims(i_sim).name).time, JRA.(sims(i_sim).name).F_res, 'LineWidth', 1.5, 'Color', 'magenta')

    hold on;
    %scatter(JRA.(sims(i_sim).name).time(max_point + round(numel(JRA.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
    xlabel('Time (s)');
    ylabel('JRF [Resultant] (N)');
    ylim(y_lims_force{4});
    hold on;


end

% Plot mean scapula data on top of other trials
figure(1)
plot(JRA.sim_M35afr903vv.time, JRA.sim_M35afr903vv.data(:, 1), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), '..\..\OpenSim\Out\Moco\Analysis\Plots\Fx.tiff','tiff')
figure(2)
plot(JRA.sim_M35afr903vv.time, JRA.sim_M35afr903vv.data(:, 2), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), '..\..\OpenSim\Out\Moco\Analysis\Plots\Fy.tiff','tiff')
figure(3)
plot(JRA.sim_M35afr903vv.time, JRA.sim_M35afr903vv.data(:, 3), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), '..\..\OpenSim\Out\Moco\Analysis\Plots\Fz.tiff','tiff')
figure(4)
plot(JRA.sim_M35afr903vv.time, JRA.sim_M35afr903vv.F_res(:, 1), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(4), '..\..\OpenSim\Out\Moco\Analysis\Plots\Fres.tiff','tiff')

close all

%% By mode
for i_mode = 1:numel(modes)

    for i_sd = 1:numel(SDs)

        morphology = [modes{i_mode} '_' SDs{i_sd}];

        idx_sim_table =  find(contains(log_table.Scapula_morphology, morphology));

        for i_joint = 1:3

            figure (i_joint)
            subplot(2,3, i_mode)
%             plot(JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).time,...
%                 JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_joint), ...
%                 'LineWidth', 1.5,...
%                 'Color', sd_colours(i_sd,:))

            plot(shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_sim_table}]),...
                JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_joint), ...
                'LineWidth', 1.5,...
                'Color', sd_colours(i_sd,:))

            hold on;
            %             scatter(JRA.(sims(i_sim).name).time(max_point + round(numel(JRA.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
            xlabel('Glenohumeral ABD (deg)');
            ylabel(['JRF [' jointF_to_plot{i_joint}(end) '] (N)']);
            ylim(y_lims_force{i_joint});
            title(modes{i_mode})
            hold on;

        end

        figure(4)
        subplot(2,3, i_mode)
%         plot(JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).time, ...
%             JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).F_res ,...
%             'LineWidth', 1.5, ...
%             'Color', sd_colours(i_sd,:))
        
        plot(shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_sim_table}]), ...
            JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).F_res ,...
            'LineWidth', 1.5, ...
            'Color', sd_colours(i_sd,:))    

        hold on;
        %scatter(JRA.(sims(i_sim).name).time(max_point + round(numel(JRA.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
        xlabel('Glenohumeral ABD (deg)','FontWeight', 'bold');
        ylabel('JRF [Resultant] (N)', 'FontWeight', 'bold');
        ylim(y_lims_force{4});
        xlim([0 80])
        title(modes{i_mode})
        hold on;

    end

    figure(1)
    plot(shoulder_elv_theta.sim_M35afr903vv, JRA.sim_M35afr903vv.data(:, 1), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
    figure(2)
    plot(shoulder_elv_theta.sim_M35afr903vv, JRA.sim_M35afr903vv.data(:, 2), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
    figure(3)
    plot(shoulder_elv_theta.sim_M35afr903vv, JRA.sim_M35afr903vv.data(:, 3), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
    figure(4)
    plot(shoulder_elv_theta.sim_M35afr903vv, JRA.sim_M35afr903vv.F_res(:, 1), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')

end

figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), '..\..\OpenSim\Out\Moco\Analysis\Plots\Fx_vs_ABD_modes.tiff','tiff')
figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), '..\..\OpenSim\Out\Moco\Analysis\Plots\Fy_vs_ABD_modes.tiff','tiff')
figure(3)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), '..\..\OpenSim\Out\Moco\Analysis\Plots\Fz_vs_ABD_modes.tiff','tiff')
figure(4)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(4), '..\..\OpenSim\Out\Moco\Analysis\Plots\Fres_vs_ABD_modes.tiff','tiff')

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
plot(TF.sim_M35afr903vv.time, TF.sim_M35afr903vv.data(:, 1), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), '..\..\OpenSim\Out\Moco\Analysis\Plots\Force_DELT1.tiff','tiff')
figure(2)
plot(TF.sim_M35afr903vv.time, TF.sim_M35afr903vv.data(:, 2), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), '..\..\OpenSim\Out\Moco\Analysis\Plots\Force_DELT2.tiff','tiff')
figure(3)
plot(TF.sim_M35afr903vv.time, TF.sim_M35afr903vv.data(:, 3), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), '..\..\OpenSim\Out\Moco\Analysis\Plots\Force_DELT3.tiff','tiff')

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
    plot(shoulder_elv_theta.sim_M35afr903vv, TF.sim_M35afr903vv.data(:, 1), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
    figure(2)
    plot(shoulder_elv_theta.sim_M35afr903vv, TF.sim_M35afr903vv.data(:, 2), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
    figure(3)
    plot(shoulder_elv_theta.sim_M35afr903vv, TF.sim_M35afr903vv.data(:, 3), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')

end


figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), '..\..\OpenSim\Out\Moco\Analysis\Plots\Force_vs_ABD_DELT1_modes.tiff','tiff')
figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), '..\..\OpenSim\Out\Moco\Analysis\Plots\Force_vs_ABD_DELT2_modes.tiff','tiff')
figure(3)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), '..\..\OpenSim\Out\Moco\Analysis\Plots\Force_vs_ABD_DELT3_modes.tiff','tiff')

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


        plot(ML.(sims(i_sim).name).time, ML.(sims(i_sim).name).data(:, i_mus), 'LineWidth', 1.5, 'Color', jrf_colours(i_mus,:))

        hold on;
        %         scatter(TF.(sims(i_sim).name).time(max_point + round(numel(TF.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
        xlabel('Time (s)');
        ylabel([' Muscle Length [' muscleF_to_plot{i_mus} '] (m)']);
        hold on;
    end

end


% Plot mean scapula data on top of other trials
figure(1)
plot(ML.sim_M35afr903vv.time, ML.sim_M35afr903vv.data(:, 1), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), '..\..\OpenSim\Out\Moco\Analysis\Plots\Length_DELT1.tiff','tiff')
figure(2)
plot(ML.sim_M35afr903vv.time, ML.sim_M35afr903vv.data(:, 2), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), '..\..\OpenSim\Out\Moco\Analysis\Plots\Length_DELT2.tiff','tiff')
figure(3)
plot(ML.sim_M35afr903vv.time, ML.sim_M35afr903vv.data(:, 3), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), '..\..\OpenSim\Out\Moco\Analysis\Plots\Length_DELT3.tiff','tiff')

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
    plot(shoulder_elv_theta.sim_M35afr903vv, ML.sim_M35afr903vv.data(:, 1), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
    figure(2)
    plot(shoulder_elv_theta.sim_M35afr903vv, ML.sim_M35afr903vv.data(:, 2), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
    figure(3)
    plot(shoulder_elv_theta.sim_M35afr903vv, ML.sim_M35afr903vv.data(:, 3), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')

end


figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), '..\..\OpenSim\Out\Moco\Analysis\Plots\Length_vs_ABD_DELT1_modes.tiff','tiff')
figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), '..\..\OpenSim\Out\Moco\Analysis\Plots\Length_vs_ABD_DELT2_modes.tiff','tiff')
figure(3)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), '..\..\OpenSim\Out\Moco\Analysis\Plots\Length_vs_ABD_DELT3_modes.tiff','tiff')
close all

%% Muscle Activations
for i_sim = 1: num_sims

    % Import ForceReporter data
    muscle_activation = importdata([analysis_folder sims(i_sim).name '/MocoSol_LateralReach.sto']);

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
plot(MAct.sim_M35afr903vv.time, MAct.sim_M35afr903vv.data(:, 1), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), '..\..\OpenSim\Out\Moco\Analysis\Plots\Activation_DELT1.tiff','tiff')
figure(2)
plot(MAct.sim_M35afr903vv.time, MAct.sim_M35afr903vv.data(:, 2), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), '..\..\OpenSim\Out\Moco\Analysis\Plots\Activation_DELT2.tiff','tiff')
figure(3)
plot(MAct.sim_M35afr903vv.time, MAct.sim_M35afr903vv.data(:, 3), 'LineWidth', 1.5, 'Color', 'black')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), '..\..\OpenSim\Out\Moco\Analysis\Plots\Activation_DELT3.tiff','tiff')

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
    plot(MAct.sim_M35afr903vv.time, MAct.sim_M35afr903vv.data(:, 1), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
    figure(2)
    plot(MAct.sim_M35afr903vv.time, MAct.sim_M35afr903vv.data(:, 2), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
    figure(3)
    plot(MAct.sim_M35afr903vv.time, MAct.sim_M35afr903vv.data(:, 3), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')

end


figure(1)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), '..\..\OpenSim\Out\Moco\Analysis\Plots\Activation_DELT1_modes.tiff','tiff')
figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), '..\..\OpenSim\Out\Moco\Analysis\Plots\Activation_DELT2_modes.tiff','tiff')
figure(3)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), '..\..\OpenSim\Out\Moco\Analysis\Plots\Activation_DELT3_modes.tiff','tiff')

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


            plot(MM.(sims(i_sim).name).time, MM.(sims(i_sim).name).data(:, i_mus), 'LineWidth', 1.5, 'Color', jrf_colours(i_mus,:))

            hold on;
            %         scatter(TF.(sims(i_sim).name).time(max_point + round(numel(TF.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
            xlabel('Time (s)');
            ylabel(['Muscle Moment [' muscleF_to_plot{i_mus} '] (Nm)']);
            hold on;
        end

    end


    % Plot mean scapula data on top of other trials
    figure(1)
    plot(MM.sim_M35afr903vv.time, MM.sim_M35afr903vv.data(:, 1), 'LineWidth', 1.5, 'Color', 'black')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(1), ['..\..\OpenSim\Out\Moco\Analysis\Plots\Moment_DELT1_' coord '.tiff'],'tiff')
    figure(2)
    plot(MM.sim_M35afr903vv.time, MM.sim_M35afr903vv.data(:, 2), 'LineWidth', 1.5, 'Color', 'black')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(2), ['..\..\OpenSim\Out\Moco\Analysis\Plots\Moment_DELT2_' coord '.tiff'],'tiff')
    figure(3)
    plot(MM.sim_M35afr903vv.time, MM.sim_M35afr903vv.data(:, 3), 'LineWidth', 1.5, 'Color', 'black')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(3), ['..\..\OpenSim\Out\Moco\Analysis\Plots\Moment_DELT3_' coord '.tiff'],'tiff')

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
        plot(shoulder_elv_theta.sim_M35afr903vv, MM.sim_M35afr903vv.data(:, 1), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
        figure(2)
        plot(shoulder_elv_theta.sim_M35afr903vv, MM.sim_M35afr903vv.data(:, 2), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
        figure(3)
        plot(shoulder_elv_theta.sim_M35afr903vv, MM.sim_M35afr903vv.data(:, 3), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')

    end


    figure(1)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(1), ['..\..\OpenSim\Out\Moco\Analysis\Plots\Moment_vs_ABD_DELT1_' coord '_modes.tiff'],'tiff')
    figure(2)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(2), ['..\..\OpenSim\Out\Moco\Analysis\Plots\Moment_vs_ABD_DELT2_' coord '_modes.tiff'],'tiff')
    figure(3)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(3), ['..\..\OpenSim\Out\Moco\Analysis\Plots\Moment_vs_ABD_DELT3_' coord '_modes.tiff'],'tiff')

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


            plot(MA.(sims(i_sim).name).time, MA.(sims(i_sim).name).data(:, i_mus), 'LineWidth', 1.5, 'Color', jrf_colours(i_mus,:))

            hold on;
            %         scatter(TF.(sims(i_sim).name).time(max_point + round(numel(TF.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
            xlabel('Time (s)');
            ylabel(['Muscle Moment Arm [' muscleF_to_plot{i_mus} '] (m)']);
            hold on;
        end

    end


    % Plot mean scapula data on top of other trials
    figure(1)
    plot(MA.sim_M35afr903vv.time, MA.sim_M35afr903vv.data(:, 1), 'LineWidth', 1.5, 'Color', 'black')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(1), ['..\..\OpenSim\Out\Moco\Analysis\Plots\MomentArm_DELT1_' coord '.tiff'],'tiff')
    figure(2)
    plot(MA.sim_M35afr903vv.time, MA.sim_M35afr903vv.data(:, 2), 'LineWidth', 1.5, 'Color', 'black')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(2), ['..\..\OpenSim\Out\Moco\Analysis\Plots\MomentArm_DELT2_' coord '.tiff'],'tiff')
    figure(3)
    plot(MA.sim_M35afr903vv.time, MA.sim_M35afr903vv.data(:, 3), 'LineWidth', 1.5, 'Color', 'black')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(3), ['..\..\OpenSim\Out\Moco\Analysis\Plots\MomentArm_DELT3_' coord '.tiff'],'tiff')

    close all

    %% By mode
    for i_mode = 1:numel(modes)

        for i_sd = 1:numel(SDs)

            morphology = [modes{i_mode} '_' SDs{i_sd}];

            idx_sim_table =  find(contains(log_table.Scapula_morphology, morphology));

            for i_mus = 1:3

                figure (i_mus)
                subplot(2,3, i_mode)
                plot(MA.(['sim_' log_table.Model_Hash{idx_sim_table}]).time,...
                    MA.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus), ...
                    'LineWidth', 1.5,...
                    'Color', sd_colours(i_sd,:))

                hold on;
                %             scatter(JRA.(sims(i_sim).name).time(max_point + round(numel(JRA.(sims(i_sim).name).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
                xlabel('Time (s)');
                ylabel(['Muscle Moment Arm [' muscleF_to_plot{i_mus} '] (m)']);
                ylim(y_lims_muscle_moment_arm{i_coord, i_mus}/1000);
                title(modes{i_mode})
                hold on;

            end


        end

        figure(1)
        plot(MA.sim_M35afr903vv.time, MA.sim_M35afr903vv.data(:, 1), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
        figure(2)
        plot(MA.sim_M35afr903vv.time, MA.sim_M35afr903vv.data(:, 2), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')
        figure(3)
        plot(MA.sim_M35afr903vv.time, MA.sim_M35afr903vv.data(:, 3), 'LineWidth', 1.5, 'Color', 'black','LineStyle','-.')

    end


    figure(1)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(1), ['..\..\OpenSim\Out\Moco\Analysis\Plots\MomentArm_DELT1_' coord '_modes.tiff'],'tiff')
    figure(2)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(2), ['..\..\OpenSim\Out\Moco\Analysis\Plots\MomentArm_DELT2_' coord '_modes.tiff'],'tiff')
    figure(3)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(3), ['..\..\OpenSim\Out\Moco\Analysis\Plots\MomentArm_DELT3_' coord '_modes.tiff'],'tiff')

    close all
end