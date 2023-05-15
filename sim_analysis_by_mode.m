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
flag_HairTouch      = false;
flag_LateralReach   = false;
flag_UpwardReach    = true;
flag_AddDummyCoord  = true;

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

axis_label_size = 18;
line_width = 2.0;

jointF_to_plot = {'unrothum_on_scapula_in_scapula_offset_fx'...
    'unrothum_on_scapula_in_scapula_offset_fy',...
    'unrothum_on_scapula_in_scapula_offset_fz'};

muscleF_to_plot = {'DELT1',...
    'DELT2',...
    'DELT3'};

coords_to_plot = {
    'shoulder_elv'
    'elv_angle'};

modes = {
%         'm4',...
        'm2',...
%         'm5',...
    %     'm6',...
    %     'm7',...
    %     'm9'
    };

% in ascending order for colour scheme
SDs = {'-3',...
    '-1',...
    '1',...
    '3'
    };

analysis_folder = '../../OpenSim/Out/Moco/Analysis/ISTA_v01/';

% List simulation folders
%sims = dir([analysis_folder '/sim_*']);

log_table = readtable([analysis_folder 'RTSA_model_log_ISTA_v01.csv']); % 'RTSA_model_log_table_final_rerun_20230404.csv'

if flag_HairTouch == true

    % Greens
    sd_colours = [229,245,224;
        161,217,155;
        65,171,93;
        0,109,44]./255;

    y_lims_force = {[-300, 300],...
        [-100, 600],...
        [0, 1500],...
        [0, 1500]};

    y_lims_q = {[0 120], [0 120];... % shoulder_elv
        [-5 20], [-40 60];... % elv_angle
        [-50 30], [-50 45]}; % shoulder_rot

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

    y_lims_muscle_moment_arm = {
        [5, 45], [20, 45], [-5, 15];
        [-10, 25], [-40, 5], [-35, 0]
        };

    task_name = 'HairTouch';

    % Get sims specific to task
    task_idx = find(contains(log_table.Task_name, task_name));
    log_table = log_table(task_idx,:);

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

    y_lims_q = {[0 90], [0 120];... % shoulder_elv
        [-40 5], [-40 60];... % elv_angle
        [-10 10], [-50 45]}; % shoulder_rot

    y_lims_muscle_force = {[0, 100],...
        [0, 400],...
        [0, 400]};

    y_lims_muscle_lengths = {[0.13, 0.22],...
        [0.10, 0.18],...
        [0.15, 0.19]};

    y_lims_muscle_moment = {%[0, 1], [-3, 0], [-2, 0];...
        [0, 5], [0, 10], [-1, 1];...
        %[-1.5, 0.2], [-7, 0], [-1.2, 0]...
        };

    y_lims_muscle_moment_arm = {
        [5, 45], [20, 55], [-5, 15];
        [0, 35], [-15, 5], [-30, 0]
        };

    task_name = 'LateralReach';

    % Get sims specific to task
    task_idx = find(contains(log_table.Task_name, task_name));
    log_table = log_table(task_idx,:);

elseif flag_UpwardReach == true

    % Blues
    sd_colours = [222, 235, 247;
        158, 202, 225;
        66, 146, 198;
        8, 69, 148]./255;

    y_lims_force = {[-400, 200],...
        [0, 500],...
        [0, 1100],...
        [0, 800]};

    y_lims_q = {[0 90], [0 120];...
        [0 70], [-40 60];...
        [0 45], [-50 45]};

    y_lims_muscle_force = {[0, 150],...
        [0, 600],...
        [0, 250]};

    y_lims_muscle_lengths = {[0.13, 0.22],...
        [0.10, 0.18],...
        [0.15, 0.19]};

    y_lims_muscle_moment = {%[0, 1], [-3, 0], [-2, 0];...
        [0, 5], [0, 10], [-1, 1];...
        %[-1.5, 0.2], [-7, 0], [-1.2, 0]...
        };

    y_lims_muscle_moment_arm = {
        [5, 45], [20, 45], [-5, 10];
        [-5, 20], [-30, 5], [-20, 0]
        };

    task_name = 'UpwardReach';

    % Get sims specific to task
    task_idx = find(contains(log_table.Task_name, task_name));
    log_table = log_table(task_idx,:);

end
num_sims = numel(log_table.Model_Hash);
idx_ave_morph = find(contains(log_table.Scapula_morphology,'m1_0'));
% sims = find(contains(sims.name,log_table.Model_Hash))

%% Re-run AnalysisTool if needed
if flag_AnalysisTool == true

    for i_sim = 1:num_sims
        model_file = ['..\..\OpenSim\In\Models\RTSA_Adjusted\FSModel_GHJoint_' log_table.Model_Hash{i_sim} '.osim'];
        osim_model = Model(model_file);

        if flag_AddDummyCoord == true
            % Add Lateral Axis coord to Shoulder
            shoulder0 = osim_model.getJointSet.get('shoulder0');
            shoulder0_dc = CustomJoint.safeDownCast(shoulder0);

            new_axis = Vec3(0,1,0);
            dummy_coord = Coordinate();
            dummy_coord.setName('dummy_coord')
            % dummy_coord_name = ArrayStr('dummy_coord');

            shoulder0_dc.set_coordinates(1,dummy_coord)

            shoulder0_dc.upd_SpatialTransform().upd_rotation3.set_coordinates(0,'dummy_coord');
            shoulder0_dc.upd_SpatialTransform().upd_rotation3.set_axis(new_axis);

            coeffs = ArrayDouble();
            coeffs.set(0,1);
            coeffs.set(1,0);

            coord_func = LinearFunction();
            % coord_func_dc = Function.safeDownCast(coord_func);
            coord_func.setCoefficients(coeffs);
            shoulder0_dc.upd_SpatialTransform().upd_rotation3.set_function(coord_func);

            % Re name file with dummy coord
            modle_file_dCoord = ['..\..\OpenSim\In\Models\RTSA_Adjusted\FSModel_GHJoint_' log_table.Model_Hash{i_sim} '_dummyCoord.osim'];

            % Re-assemble model
            osim_model.finalizeConnections()
            osim_model.print(modle_file_dCoord);

        end


        solution_file = [analysis_folder ['sim_' log_table.Model_Hash{i_sim}] '/MocoSol_' task_name '.sto'];

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
        analyzeTool=AnalyzeTool('..\..\OpenSim\In\Setup_files\Analysis\template_MA_dummyCoord.xml',0);
        analyzeTool.setName('Moco');
        analyzeTool.setInitialTime(0);
        analyzeTool.setFinalTime(time_array.getLast);
        %             analyzeTool.setStatesStorage(states_storage);
        analyzeTool.setStatesFileName(solution_file);
        %             analyzeTool.setModel(osim_model);
        if flag_AddDummyCoord == true
            analyzeTool.setModelFilename(modle_file_dCoord)
        else
            analyzeTool.setModelFilename(model_file)
        end

        analyzeTool.setResultsDir(['..\..\OpenSim\Out\Moco\Analysis\paper_v02\Final\' ['sim_' log_table.Model_Hash{i_sim}] '\']);

        %         JR_Analysis = analyzeTool.updAnalysisSet.get(0);
        %         JR_downCast = JointReaction.safeDownCast(JR_Analysis);
        %         JR_downCast.setEndTime(time_array.getLast);
        %         JR_downCast.setInFrame(in_frames);
        %         JR_downCast.setOnBody(on_bodies);
        %         JR_downCast.setJointNames(joints);
        %
        %         FR = analyzeTool.updAnalysisSet.get(1);
        %         FR.setEndTime(time_array.getLast);
        %
        %         MA = analyzeTool.updAnalysisSet.get(2);
        %         MA.setEndTime(time_array.getLast);

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
    joint_reaction = importdata([analysis_folder ['sim_' log_table.Model_Hash{i_sim}] '/Moco_JointReaction_ReactionLoads.sto']);
    % Import States and Kinematics
    states = importdata([analysis_folder ['sim_' log_table.Model_Hash{i_sim}] '/MocoSol_' task_name '.sto']);
    shoulder_elv_theta.(['sim_' log_table.Model_Hash{i_sim}]) = rad2deg(states.data(:,19));
    elv_angle_theta.(['sim_' log_table.Model_Hash{i_sim}]) = rad2deg(states.data(:,18));
    shoulder_rot_theta.(['sim_' log_table.Model_Hash{i_sim}]) = rad2deg(states.data(:,21));

    % Import Kinematics data
    %%% sim_JRF = importdata([analysis_folder ['sim_' log_table.Model_Hash{i_sim}] '/MocoSol_' task_name '.sto']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% JRF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get time
    if flag_NormTime == true
        time_vec = linspace(0, 100, 101)';
        JRA.(['sim_' log_table.Model_Hash{i_sim}]).time_og = 100*(joint_reaction.data(:,1)/max(joint_reaction.data(:,1)));
        JRA.(['sim_' log_table.Model_Hash{i_sim}]).time = 100*(joint_reaction.data(:,1)/max(joint_reaction.data(:,1)));
    else
        JRA.(['sim_' log_table.Model_Hash{i_sim}]).time = joint_reaction.data(:,1);
    end


    for i_joint = 1:numel(jointF_to_plot)
        % Get handle of coordinates to plot
        JRA.(['sim_' log_table.Model_Hash{i_sim}]).pos(:, i_joint) = find(contains(joint_reaction.colheaders, jointF_to_plot{i_joint}));
        JRA.(['sim_' log_table.Model_Hash{i_sim}]).label{1,i_joint} = joint_reaction.colheaders{JRA.(['sim_' log_table.Model_Hash{i_sim}]).pos(i_joint)};

        if flag_NormTime == false
            if contains(JRA.(['sim_' log_table.Model_Hash{i_sim}]).label{1,i_joint}, 'fz')
                JRA.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_joint) = abs(joint_reaction.data(:,JRA.(['sim_' log_table.Model_Hash{i_sim}]).pos(i_joint)));
            else
                JRA.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_joint) = joint_reaction.data(:,JRA.(['sim_' log_table.Model_Hash{i_sim}]).pos(i_joint));
            end
        elseif flag_NormTime == true

            if contains(JRA.(['sim_' log_table.Model_Hash{i_sim}]).label{1,i_joint}, 'fz')
                JRA.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_joint) = interp1(JRA.(['sim_' log_table.Model_Hash{i_sim}]).time_og,...
                    abs(joint_reaction.data(:,JRA.(['sim_' log_table.Model_Hash{i_sim}]).pos(i_joint))),...
                    time_vec);
            else
                JRA.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_joint) = interp1(JRA.(['sim_' log_table.Model_Hash{i_sim}]).time_og,...
                    joint_reaction.data(:,JRA.(['sim_' log_table.Model_Hash{i_sim}]).pos(i_joint)),...
                    time_vec);
            end

            JRA.(['sim_' log_table.Model_Hash{i_sim}]).time = time_vec;


        end

        % Create JRF plots
        figure (i_joint)

        plot(JRA.(['sim_' log_table.Model_Hash{i_sim}]).time, JRA.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_joint), 'LineWidth', line_width, 'Color', jrf_colours(i_joint,:))

        hold on;
        xlabel('Movement Duration (%)', 'FontWeight', 'bold', 'FontSize', axis_label_size);
        ylabel(['JRF [' jointF_to_plot{i_joint}(end) '] (N)'], 'FontWeight', 'bold', 'FontSize', axis_label_size);
        ylim(y_lims_force{i_joint});
        hold on;



    end

    % Normalise Kinematics
    shoulder_elv_theta.(['sim_' log_table.Model_Hash{i_sim}]) = interp1(JRA.(['sim_' log_table.Model_Hash{i_sim}]).time_og,...
        shoulder_elv_theta.(['sim_' log_table.Model_Hash{i_sim}]),...
        time_vec);

    elv_angle_theta.(['sim_' log_table.Model_Hash{i_sim}]) = interp1(JRA.(['sim_' log_table.Model_Hash{i_sim}]).time_og,...
        elv_angle_theta.(['sim_' log_table.Model_Hash{i_sim}]),...
        time_vec);

    shoulder_rot_theta.(['sim_' log_table.Model_Hash{i_sim}]) = interp1(JRA.(['sim_' log_table.Model_Hash{i_sim}]).time_og,...
        shoulder_rot_theta.(['sim_' log_table.Model_Hash{i_sim}]),...
        time_vec);

    % Calculate resultant Force
    JRA.(['sim_' log_table.Model_Hash{i_sim}]).F_res = sqrt(JRA.(['sim_' log_table.Model_Hash{i_sim}]).data(:, 1).^2 + JRA.(['sim_' log_table.Model_Hash{i_sim}]).data(:, 2).^2 + JRA.(['sim_' log_table.Model_Hash{i_sim}]).data(:, 3).^2);

    figure(4)

    plot(JRA.(['sim_' log_table.Model_Hash{i_sim}]).time, JRA.(['sim_' log_table.Model_Hash{i_sim}]).F_res, 'LineWidth', line_width, 'Color', jrf_colours(4,:))

    hold on;
    %scatter(JRA.(['sim_' log_table.Model_Hash{i_sim}]).time(max_point + round(numel(JRA.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
    xlabel('Movement Duration (%)','FontWeight', 'bold','FontSize', axis_label_size);
    ylabel('Resultant JRF (N)','FontWeight', 'bold','FontSize', axis_label_size);
    ylim(y_lims_force{4});
    hold on;

    % Plot Joint Kinematics

    figure (11)

    plot(JRA.(['sim_' log_table.Model_Hash{i_sim}]).time, shoulder_elv_theta.(['sim_' log_table.Model_Hash{i_sim}]), 'LineWidth', line_width, 'Color', jrf_colours(1,:))

    hold on;
    ylim(y_lims_q{1,2});
    xlabel('Movement Duration (%)','FontWeight', 'bold','FontSize', axis_label_size);
    ylabel(['Shoulder Elevation (' char(176) ')'],'FontWeight', 'bold','FontSize', axis_label_size);
    hold on;

    figure (12)

    plot(JRA.(['sim_' log_table.Model_Hash{i_sim}]).time, elv_angle_theta.(['sim_' log_table.Model_Hash{i_sim}]), 'LineWidth', line_width, 'Color', jrf_colours(2,:))

    hold on;
    ylim(y_lims_q{2,2});
    xlabel('Movement Duration (%)','FontWeight', 'bold','FontSize', axis_label_size);
    ylabel(['Elevation Plane (' char(176) ')'],'FontWeight', 'bold','FontSize', axis_label_size);
    hold on;

    figure (13)

    plot(JRA.(['sim_' log_table.Model_Hash{i_sim}]).time, shoulder_rot_theta.(['sim_' log_table.Model_Hash{i_sim}]), 'LineWidth', line_width, 'Color', jrf_colours(3,:))

    hold on;
    ylim(y_lims_q{3,2});
    xlabel('Movement Duration (%)','FontWeight', 'bold','FontSize', axis_label_size);
    ylabel(['Shoulder Rotation (' char(176) ')'],'FontWeight', 'bold','FontSize', axis_label_size);
    hold on;

end

% Plot mean scapula data on top of other trials
figure(1)
set(gca,'FontWeight', 'bold')
plot(JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 1), 'LineWidth', line_width, 'Color', 'black')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), [analysis_folder task_name '/Fx_norm.tiff'],'tiff')

figure(2)
set(gca,'FontWeight', 'bold')
plot(JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 2), 'LineWidth', line_width, 'Color', 'black')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), [analysis_folder task_name '/Fy_norm.tiff'],'tiff')

figure(3)
set(gca,'FontWeight', 'bold')
plot(JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 3), 'LineWidth', line_width, 'Color', 'black')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), [analysis_folder task_name '/Fz_norm.tiff'],'tiff')

figure(4)
set(gca,'FontWeight', 'bold')
plot(JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).F_res(:, 1), 'LineWidth', line_width, 'Color', 'black')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(4), [analysis_folder task_name '/Fres_norm.tiff'],'tiff')

figure(11)
set(gca,'FontWeight', 'bold')
plot(JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_ave_morph}]), 'LineWidth', line_width, 'Color', 'black')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(11), [analysis_folder task_name '/shoulder_elv_norm.tiff'],'tiff')

figure(12)
set(gca,'FontWeight', 'bold')
plot(JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, elv_angle_theta.(['sim_' log_table.Model_Hash{idx_ave_morph}]), 'LineWidth', line_width, 'Color', 'black')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);['sim_' log_table.Model_Hash{idx_sim_table}]
saveas(figure(12), [analysis_folder task_name '/elv_ang_norm.tiff'],'tiff')

figure(13)
set(gca,'FontWeight', 'bold')
plot(JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, shoulder_rot_theta.(['sim_' log_table.Model_Hash{idx_ave_morph}]), 'LineWidth', line_width, 'Color', 'black')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(13), [analysis_folder task_name '/shoulder_rot_norm.tiff'],'tiff')

save([analysis_folder task_name '/data_JRA'], 'JRA')
save([analysis_folder task_name '/data_shoulder_elv_theta'], 'shoulder_elv_theta')
save([analysis_folder task_name '/data_elv_angle_theta'], "elv_angle_theta")
save([analysis_folder task_name '/data_shoulder_rot_theta'], "shoulder_rot_theta")
close all

%% By mode
for i_mode = 1:numel(modes)

    for i_sd = 1:numel(SDs)

        morphology = [modes{i_mode} '_' SDs{i_sd}];

        % If comparing two methods of the same mode (i.e. different trials)
        n_idx_sim_table =  find(contains(log_table.Scapula_morphology, morphology));

        for i_trial = 1:numel(n_idx_sim_table)

            idx_sim_table = n_idx_sim_table(i_trial);

            if log_table.Date_Generated(idx_sim_table)>=datetime(2023,5,1) && log_table.Date_Generated(idx_sim_table)<=datetime(2023,5,10)
                line_style = '-';
                % Reds
                sd_colours = [254,153,41;
                    236,112,20;
                    204,76,2;
                    140,45,4]./255;

                F_res_plus5(:,i_sd) = JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).F_res;
            
            elseif log_table.Date_Generated(idx_sim_table)>=datetime(2023,5,10)
                line_style = '-';
                % Pinks
                sd_colours = [221,52,151;
                    174,1,126;
                    122,1,119;
                    73,0,106]./255;

                F_res_plus10(:,i_sd) = JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).F_res;
            else
                line_style = '-';
                % Blues
                sd_colours = [222, 235, 247;
                    158, 202, 225;
                    66, 146, 198;
                    8, 69, 148]./255;

                F_res(:,i_sd) = JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).F_res;
            end

            for i_joint = 1:3

                figure (i_joint)
                %             subplot(2,3, i_mode) % subplot(6,3, i_mode + (i_mode -1)*2)
                plot(JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).time,...
                    JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_joint), ...
                    'LineWidth', line_width,...
                    'Color', sd_colours(i_sd,:),...
                    'LineStyle', line_style);

                %             plot(shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_sim_table}]),...
                %                 JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_joint), ...
                %                 'LineWidth', line_width,...
                %                 'Color', sd_colours(i_sd,:))

                hold on;
                %             scatter(JRA.(['sim_' log_table.Model_Hash{i_sim}]).time(max_point + round(numel(JRA.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
                xlabel('Movement Duration (%)', 'FontWeight', 'bold','FontSize', axis_label_size);
                ylabel(['JRF [' jointF_to_plot{i_joint}(end) '] (N)'], 'FontWeight', 'bold','FontSize', axis_label_size);
                ylim(y_lims_force{i_joint});
                title(modes{i_mode})
                hold on;

            end

            figure(4)
            %         subplot(2,3, i_mode)
            plot(JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).time, ...
                JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).F_res ,...
                'LineWidth', line_width, ...
                'Color', sd_colours(i_sd,:),...
                'LineStyle', line_style);

            hold on;
            %scatter(JRA.(['sim_' log_table.Model_Hash{i_sim}]).time(max_point + round(numel(JRA.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
            xlabel('Movement Duration (%)','FontWeight', 'bold','FontSize', axis_label_size);
            ylabel('Resultant JRF (N)', 'FontWeight', 'bold','FontSize', axis_label_size);
            ylim(y_lims_force{4});
            %xlim([0 80])
            title(modes{i_mode})
            hold on;

            % Plot Force ratio
            figure(5)

            fr = sqrt(JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, 3).^2)./...
                (sqrt(JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, 1).^2) + sqrt(JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, 2).^2));


            %         subplot(2,3, i_mode)
            plot(JRA.(['sim_' log_table.Model_Hash{i_sim}]).time,...
                fr,...
                'LineWidth', line_width, ...
                'Color', sd_colours(i_sd,:), ...
                'LineStyle', line_style)
            hold on
            ylabel('FR', 'FontWeight', 'bold','FontSize', axis_label_size);
            xlabel('Movement Duration (%)','FontWeight', 'bold','FontSize', axis_label_size);
            ylim([0, 1.6])
            title(modes{i_mode})


            % Plot shoulder_elv by mode
            figure(6)
            %         subplot(2,3, i_mode)

            plot(JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).time, ...
                shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_sim_table}]), ...
                'LineWidth', line_width, ...
                'Color', sd_colours(i_sd,:),...
                'LineStyle', line_style)

            hold on;
            xlabel('Movement Duration (%)','FontWeight', 'bold','FontSize', axis_label_size);
            ylabel(['Shoulder Elevation (' char(176) ')'], 'FontWeight', 'bold','FontSize', axis_label_size);
            ylim(y_lims_q{1,1});
            title(modes{i_mode})
            hold on;

            % Plot elv_angle by mode
            figure(7)
            %         subplot(2,3, i_mode)

            plot(JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).time, ...
                elv_angle_theta.(['sim_' log_table.Model_Hash{idx_sim_table}]), ...
                'LineWidth', line_width, ...
                'Color', sd_colours(i_sd,:), ...
                'LineStyle', line_style)

            hold on;
            xlabel('Movement Duration (%)','FontWeight', 'bold','FontSize', axis_label_size);
            ylabel(['Elevation Plane (' char(176) ')'], 'FontWeight', 'bold','FontSize', axis_label_size);
            ylim(y_lims_q{2,1});
            title(modes{i_mode})
            hold on;

            % Plot shoulder_rot by mode
            figure(8)
            %         subplot(2,3, i_mode)

            plot(JRA.(['sim_' log_table.Model_Hash{idx_sim_table}]).time, ...
                shoulder_rot_theta.(['sim_' log_table.Model_Hash{idx_sim_table}]), ...
                'LineWidth', line_width, ...
                'Color', sd_colours(i_sd,:), ...
                'LineStyle', line_style)

            hold on;
            xlabel('Movement Duration (%)','FontWeight', 'bold','FontSize', axis_label_size);
            ylabel(['Shoulder Rotation (' char(176) ')'], 'FontWeight', 'bold','FontSize', axis_label_size);
            ylim(y_lims_q{3,1});
            title(modes{i_mode})
            hold on;

        end
    end

    %% SPM 
    % inference:
    alpha      = 0.05;
    nTests     = 3;
    p_critical = spm1d.util.p_critical_bonf(alpha, nTests);

    figure(33)
    subplot(3,1,1)
    t_baseV5 = spm1d.stats.ttest_paired(F_res',F_res_plus5');
    t_baseV5i = t_baseV5.inference(p_critical,'two_tailed', true, 'interp', true);
    t_baseV5i.plot()
    ylim([-100,25]);
    xticks([])
    title('Baseline vs 5 mm','FontSize', 12)
%     xlabel('Movement Duration (%)','FontWeight', 'bold','FontSize', axis_label_size);
    

    figure(33)
    subplot(3,1,2)
    t_baseV10 = spm1d.stats.ttest_paired(F_res',F_res_plus10');
    baseV10i = t_baseV10.inference(p_critical,'two_tailed', true, 'interp', true);
    baseV10i.plot()
    ylim([-100,25]);
    xticks([])
    title('Baseline vs 10 mm','FontSize', 12)
%     xlabel('Movement Duration (%)','FontWeight', 'bold','FontSize', axis_label_size);

    figure(33)
    subplot(3,1,3)
    t_5V10 = spm1d.stats.ttest_paired(F_res_plus5',F_res_plus10');
    t_5V10i = t_5V10.inference(p_critical,'two_tailed', true, 'interp', true);
    t_5V10i.plot()
    ylim([-100,25]);
    title('5 mm vs 10 mm','FontSize', 12)
    xlabel('Movement Duration (%)','FontWeight', 'bold','FontSize', axis_label_size);


    saveas(figure(33), [analysis_folder task_name '/SPM_mode_' modes{i_mode} '.tiff'],'tiff')

    A = [F_res'; F_res_plus5'; F_res_plus10'];
    SUB = repmat([1:4]',3,1);
    condition = [0,0,0,0,1,1,1,1,2,2,2,2]';
    figure(35)
    F  = spm1d.stats.anova1(A,condition, SUB);
    Fi = F.inference(0.05);
    Fi.plot()
    xlabel('Movement Duration (%)','FontWeight', 'bold','FontSize', axis_label_size);


    % Average morphology data
    figure(1)
    plot(JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 1), 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
    figure(2)
    plot(JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 2), 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
    figure(3)
    plot(JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 3), 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
    figure(4)
    plot(JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).F_res(:, 1), 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')

    figure(5)
    fr_av = sqrt(JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 3).^2)./...
        (sqrt(JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 1).^2) + sqrt(JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 2).^2));
    plot(JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, fr_av, 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')

    figure(6)
    plot(JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_ave_morph}]), 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
    figure(7)
    plot(JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, elv_angle_theta.(['sim_' log_table.Model_Hash{idx_ave_morph}]), 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
    figure(8)
    plot(JRA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, shoulder_rot_theta.(['sim_' log_table.Model_Hash{idx_ave_morph}]), 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')


    figure(1)
    set(gca,'FontWeight', 'bold')
    %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(1), [analysis_folder task_name '/Fx_vs_NORM_mode_' modes{i_mode} '.tiff'],'tiff')
    figure(2)
    set(gca,'FontWeight', 'bold')
    %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(2), [analysis_folder task_name '/Fy_vs_NORM_mode_' modes{i_mode} '.tiff'],'tiff')
    figure(3)
    set(gca,'FontWeight', 'bold')
    %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(3), [analysis_folder task_name '/Fz_vs_NORM_mode_' modes{i_mode} '.tiff'],'tiff')
    figure(4)
    set(gca,'FontWeight', 'bold')
    %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(4), [analysis_folder task_name '/Fres_vs_NORM_mode_' modes{i_mode} '.tiff'],'tiff')
    figure(5)
    set(gca,'FontWeight', 'bold')
    %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(5), [analysis_folder task_name '/FR_vs_NORM_mode_' modes{i_mode} '.tiff'],'tiff')
    figure(6)
    set(gca,'FontWeight', 'bold')
    %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(6), [analysis_folder task_name '/shoulder_elv_norm_mode_' modes{i_mode} '.tiff'],'tiff')
    figure(7)
    set(gca,'FontWeight', 'bold')
    %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(7), [analysis_folder task_name '/elv_angle_norm_mode_' modes{i_mode} '.tiff'],'tiff')
    figure(8)
    set(gca,'FontWeight', 'bold')
    %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(8), [analysis_folder task_name '/shoulder_rotation_norm_mode_' modes{i_mode} '.tiff'],'tiff')

    close all

end


%% Analysis muscle forces
for i_sim = 1: num_sims

    % Import ForceReporter data
    tendon_force = importdata([analysis_folder ['sim_' log_table.Model_Hash{i_sim}] '/Moco_MuscleAnalysis_FiberForce.sto']);

    % Get time
    TF.(['sim_' log_table.Model_Hash{i_sim}]).time = tendon_force.data(:,1);


    for i_mus = 1:numel(muscleF_to_plot)
        % Get handle of coordinates to plot
        TF.(['sim_' log_table.Model_Hash{i_sim}]).pos(:, i_mus) = find(contains(tendon_force.colheaders, muscleF_to_plot{i_mus}));
        TF.(['sim_' log_table.Model_Hash{i_sim}]).label{1,i_mus} = tendon_force.colheaders{TF.(['sim_' log_table.Model_Hash{i_sim}]).pos(i_mus)};


        %%%%%%%%%%%%%%% To normalise %%%%%%%%%%%%%%%%%%%%%5
        if flag_NormTime == true
            TF.(['sim_' log_table.Model_Hash{i_sim}]).time_og = 100*(tendon_force.data(:,1)/max(tendon_force.data(:,1)));
            TF.(['sim_' log_table.Model_Hash{i_sim}]).time = 100*(tendon_force.data(:,1)/max(tendon_force.data(:,1)));

            TF.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_mus) = interp1(TF.(['sim_' log_table.Model_Hash{i_sim}]).time_og,...
                tendon_force.data(:,TF.(['sim_' log_table.Model_Hash{i_sim}]).pos(i_mus)),...
                time_vec);
        else
            TF.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_mus) = tendon_force.data(:,TF.(['sim_' log_table.Model_Hash{i_sim}]).pos(i_mus));

        end

        TF.(['sim_' log_table.Model_Hash{i_sim}]).time = time_vec;

        % Create q plots
        figure (i_mus)


        plot(TF.(['sim_' log_table.Model_Hash{i_sim}]).time, TF.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_mus), 'LineWidth', line_width, 'Color', jrf_colours(i_mus,:))

        hold on;
        %         scatter(TF.(['sim_' log_table.Model_Hash{i_sim}]).time(max_point + round(numel(TF.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
        xlabel('Movement Duration (%)', 'FontWeight', 'bold','FontSize', axis_label_size);
        ylabel('Tendon Force (N)', 'FontWeight', 'bold','FontSize', axis_label_size);
        hold on;
    end

end


% Plot mean scapula data on top of other trials
figure(1)
set(gca,'FontWeight', 'bold')
plot(TF.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, TF.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 1), 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), [analysis_folder task_name '/Force_DELT1.tiff'],'tiff')
figure(2)
set(gca,'FontWeight', 'bold')
plot(TF.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, TF.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 2), 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), [analysis_folder task_name '/Force_DELT2.tiff'],'tiff')
figure(3)
set(gca,'FontWeight', 'bold')
plot(TF.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, TF.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 3), 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), [analysis_folder task_name '/Force_DELT3.tiff'],'tiff')

save([analysis_folder task_name '/data_TF'], 'TF')
close all

%% By mode
for i_mus = 1:3
    for i_mode = 1:numel(modes)
        for i_sd = 1:numel(SDs)

            morphology = [modes{i_mode} '_' SDs{i_sd}];


            % If comparing two methods of the same mode (i.e. different trials)
            n_idx_sim_table =  find(contains(log_table.Scapula_morphology, morphology));

            for i_trial = 1:numel(n_idx_sim_table)

                idx_sim_table = n_idx_sim_table(i_trial);

                if log_table.Date_Generated(idx_sim_table)>=datetime(2023,5,1) && log_table.Date_Generated(idx_sim_table)<=datetime(2023,5,10)
                    line_style = '-';
                    % Reds
                    sd_colours = [254,153,41;
                        236,112,20;
                        204,76,2;
                        140,45,4]./255;
                elseif log_table.Date_Generated(idx_sim_table)>=datetime(2023,5,10)
                    line_style = '-';
                    % Pinks
                    sd_colours = [221,52,151;
                        174,1,126;
                        122,1,119;
                        73,0,106]./255;
                else
                    line_style = '-';
                    % Blues
                    sd_colours = [222, 235, 247;
                        158, 202, 225;
                        66, 146, 198;
                        8, 69, 148]./255;
                end

                figure(1)
                plot(TF.(['sim_' log_table.Model_Hash{i_sim}]).time,...
                    TF.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus), ...
                    'LineWidth', line_width,...
                    'Color', sd_colours(i_sd,:), ...
                    'LineStyle', line_style)

                hold on;
            end
        end

        figure(1)
        plot(TF.(['sim_' log_table.Model_Hash{i_sim}]).time, TF.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, i_mus), 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')

        xlabel('Movement Duration (%)','FontWeight', 'bold', 'FontSize', axis_label_size);
        ylabel('Tendon Force (N)', 'FontWeight', 'bold', 'FontSize', axis_label_size);
        ylim(y_lims_muscle_force{i_mus});
        title(modes{i_mode})
        hold on;

        figure(1)
        set(gca,'FontWeight', 'bold')
        %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        saveas(figure(1), [analysis_folder task_name '/Force_vs_ABD_DELT' num2str(i_mus) '_' num2str(modes{i_mode}) '.tiff'],'tiff')
        close all
    end


end


% % % % % % % % % %% Muscle Lengths
% % % % % % % % % for i_sim = 1: num_sims
% % % % % % % % %
% % % % % % % % %     % Import ForceReporter data
% % % % % % % % %     muscle_length = importdata([analysis_folder ['sim_' log_table.Model_Hash{i_sim}] '/Moco_MuscleAnalysis_Length.sto']);
% % % % % % % % %
% % % % % % % % %     % Get time
% % % % % % % % %     ML.(['sim_' log_table.Model_Hash{i_sim}]).time = muscle_length.data(:,1);
% % % % % % % % %
% % % % % % % % %     for i_mus = 1:numel(muscleF_to_plot)
% % % % % % % % %         % Get handle of coordinates to plot
% % % % % % % % %         ML.(['sim_' log_table.Model_Hash{i_sim}]).pos(:, i_mus) = find(contains(muscle_length.colheaders, muscleF_to_plot{i_mus}));
% % % % % % % % %         ML.(['sim_' log_table.Model_Hash{i_sim}]).label{1,i_mus} = muscle_length.colheaders{ML.(['sim_' log_table.Model_Hash{i_sim}]).pos(i_mus)};
% % % % % % % % %
% % % % % % % % %         ML.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_mus) = muscle_length.data(:,ML.(['sim_' log_table.Model_Hash{i_sim}]).pos(i_mus));
% % % % % % % % %
% % % % % % % % %
% % % % % % % % %         % Create q plots
% % % % % % % % %         figure (i_mus)
% % % % % % % % %
% % % % % % % % %
% % % % % % % % %         plot(shoulder_elv_theta.(['sim_' log_table.Model_Hash{i_sim}]), ML.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_mus), 'LineWidth', line_width, 'Color', jrf_colours(i_mus,:))
% % % % % % % % %
% % % % % % % % %         hold on;
% % % % % % % % %         %         scatter(TF.(['sim_' log_table.Model_Hash{i_sim}]).time(max_point + round(numel(TF.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
% % % % % % % % %         xlabel('Time (s)');
% % % % % % % % %         ylabel([' Muscle Length [' muscleF_to_plot{i_mus} '] (m)']);
% % % % % % % % %         hold on;
% % % % % % % % %     end
% % % % % % % % %
% % % % % % % % % end
% % % % % % % % %
% % % % % % % % %
% % % % % % % % % % Plot mean scapula data on top of other trials
% % % % % % % % % figure(1)
% % % % % % % % % plot(ML.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, ML.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 1), 'LineWidth', line_width, 'Color', 'black')
% % % % % % % % % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % % % % % % % % saveas(figure(1), [analysis_folder task_name '/Length_DELT1.tiff'],'tiff')
% % % % % % % % % figure(2)
% % % % % % % % % plot(ML.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, ML.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 2), 'LineWidth', line_width, 'Color', 'black')
% % % % % % % % % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % % % % % % % % saveas(figure(2), [analysis_folder task_name '/Length_DELT2.tiff'],'tiff')
% % % % % % % % % figure(3)
% % % % % % % % % plot(ML.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, ML.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 3), 'LineWidth', line_width, 'Color', 'black')
% % % % % % % % % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % % % % % % % % saveas(figure(3), [analysis_folder task_name '/Length_DELT3.tiff'],'tiff')
% % % % % % % % %
% % % % % % % % % close all
% % % % % % % % %
% % % % % % % % % %% By mode
% % % % % % % % %
% % % % % % % % % for i_mode = 1:numel(modes)
% % % % % % % % %
% % % % % % % % %     for i_sd = 1:numel(SDs)
% % % % % % % % %
% % % % % % % % %         morphology = [modes{i_mode} '_' SDs{i_sd}];
% % % % % % % % %
% % % % % % % % %         idx_sim_table =  find(contains(log_table.Scapula_morphology, morphology));
% % % % % % % % %
% % % % % % % % %         for i_mus = 1:3
% % % % % % % % %
% % % % % % % % %             figure (i_mus)
% % % % % % % % %             subplot(2,3, i_mode)
% % % % % % % % %             %             plot(ML.(['sim_' log_table.Model_Hash{idx_sim_table}]).time,...
% % % % % % % % %             %                 ML.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus), ...
% % % % % % % % %             %                 'LineWidth', line_width,...
% % % % % % % % %             %                 'Color', sd_colours(i_sd,:))
% % % % % % % % %
% % % % % % % % %             plot(shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_sim_table}]),...
% % % % % % % % %                 ML.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus), ...
% % % % % % % % %                 'LineWidth', line_width,...
% % % % % % % % %                 'Color', sd_colours(i_sd,:))
% % % % % % % % %
% % % % % % % % %             hold on;
% % % % % % % % %             %             scatter(JRA.(['sim_' log_table.Model_Hash{i_sim}]).time(max_point + round(numel(JRA.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
% % % % % % % % %             xlabel(['Glenohumeral ABD (' char(176) ')']);
% % % % % % % % %             ylabel(['Muscle Length [' muscleF_to_plot{i_mus} '] (m)']);
% % % % % % % % %             ylim(y_lims_muscle_lengths{i_mus});
% % % % % % % % %             title(modes{i_mode})
% % % % % % % % %             hold on;
% % % % % % % % %
% % % % % % % % %         end
% % % % % % % % %
% % % % % % % % %
% % % % % % % % %     end
% % % % % % % % %
% % % % % % % % %     figure(1)
% % % % % % % % %     plot(shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_ave_morph}]), ML.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 1), 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
% % % % % % % % %     figure(2)
% % % % % % % % %     plot(shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_ave_morph}]), ML.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 2), 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
% % % % % % % % %     figure(3)
% % % % % % % % %     plot(shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_ave_morph}]), ML.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 3), 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
% % % % % % % % %
% % % % % % % % % end
% % % % % % % % %
% % % % % % % % %
% % % % % % % % % figure(1)
% % % % % % % % % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % % % % % % % % saveas(figure(1), [analysis_folder task_name '/Length_vs_ABD_DELT1_modes.tiff'],'tiff')
% % % % % % % % % figure(2)
% % % % % % % % % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % % % % % % % % saveas(figure(2), [analysis_folder task_name '/Length_vs_ABD_DELT2_modes.tiff'],'tiff')
% % % % % % % % % figure(3)
% % % % % % % % % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % % % % % % % % saveas(figure(3), [analysis_folder task_name '/Length_vs_ABD_DELT3_modes.tiff'],'tiff')
% % % % % % % % % close all
% % % % % % % % %
%% Muscle Activations
for i_sim = 1: num_sims

    % Import ForceReporter data
    muscle_activation = importdata([analysis_folder ['sim_' log_table.Model_Hash{i_sim}] '/MocoSol_' task_name '.sto']);

    % Get time
    MAct.(['sim_' log_table.Model_Hash{i_sim}]).time = muscle_activation.data(:,1);

    for i_mus = 1:numel(muscleF_to_plot)
        % Get handle of coordinates to plot
        MAct.(['sim_' log_table.Model_Hash{i_sim}]).pos(:, i_mus) = find(contains(muscle_activation.colheaders, [muscleF_to_plot{i_mus} '/activation']));
        MAct.(['sim_' log_table.Model_Hash{i_sim}]).label{1,i_mus} = muscle_activation.colheaders{MAct.(['sim_' log_table.Model_Hash{i_sim}]).pos(i_mus)};

        %%%%%%%%%%%%%%% To normalise %%%%%%%%%%%%%%%%%%%%%5
        if flag_NormTime == true
            MAct.(['sim_' log_table.Model_Hash{i_sim}]).time_og = 100*(muscle_activation.data(:,1)/max(muscle_activation.data(:,1)));
            MMActA.(['sim_' log_table.Model_Hash{i_sim}]).time = 100*(muscle_activation.data(:,1)/max(muscle_activation.data(:,1)));

            MAct.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_mus) = interp1(MAct.(['sim_' log_table.Model_Hash{i_sim}]).time_og,...
                muscle_activation.data(:,MAct.(['sim_' log_table.Model_Hash{i_sim}]).pos(i_mus)),...
                time_vec);
        else
            MAct.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_mus) = muscle_activation.data(:,MAct.(['sim_' log_table.Model_Hash{i_sim}]).pos(i_mus));

        end

        MAct.(['sim_' log_table.Model_Hash{i_sim}]).time = time_vec;

        % Create q plots
        figure (i_mus)


        plot(MAct.(['sim_' log_table.Model_Hash{i_sim}]).time, MAct.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_mus) *100, 'LineWidth', line_width, 'Color', jrf_colours(i_mus,:))

        hold on;
        %         scatter(TF.(['sim_' log_table.Model_Hash{i_sim}]).time(max_point + round(numel(TF.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
        xlabel('Movement Duration (%)','FontWeight', 'bold', 'FontSize', axis_label_size);
        ylabel('Muscle Activation (%)' , 'FontWeight', 'bold', 'FontSize', axis_label_size);
        ylim([0, 60]);
        hold on;
    end

end


% Plot mean scapula data on top of other trials
figure(1)
set(gca,'FontWeight', 'bold')
plot(MAct.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, MAct.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 1)*100, 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(1), [analysis_folder task_name '/Activation_DELT1.tiff'],'tiff')
figure(2)
set(gca,'FontWeight', 'bold')
plot(MAct.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, MAct.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 2)*100, 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(2), [analysis_folder task_name '/Activation_DELT2.tiff'],'tiff')
figure(3)
set(gca,'FontWeight', 'bold')
plot(MAct.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, MAct.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 3)*100, 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(figure(3), [analysis_folder task_name '/Activation_DELT3.tiff'],'tiff')

save([analysis_folder task_name '/data_MAct'], 'MAct')
close all

%% By mode
for i_mus = 1:3
    for i_mode = 1:numel(modes)

        for i_sd = 1:numel(SDs)

            morphology = [modes{i_mode} '_' SDs{i_sd}];

            % If comparing two methods of the same mode (i.e. different trials)
            n_idx_sim_table =  find(contains(log_table.Scapula_morphology, morphology));

            for i_trial = 1:numel(n_idx_sim_table)

                idx_sim_table = n_idx_sim_table(i_trial);

                if log_table.Date_Generated(idx_sim_table)>=datetime(2023,5,1) && log_table.Date_Generated(idx_sim_table)<=datetime(2023,5,10)
                    line_style = '-';
                    % Reds
                    sd_colours = [254,153,41;
                        236,112,20;
                        204,76,2;
                        140,45,4]./255;
                elseif log_table.Date_Generated(idx_sim_table)>=datetime(2023,5,10)
                line_style = '-';
                % Pinks
                sd_colours = [221,52,151;
                    174,1,126;
                    122,1,119;
                    73,0,106]./255;
                else
                    line_style = '-';
                    % Blues
                    sd_colours = [222, 235, 247;
                        158, 202, 225;
                        66, 146, 198;
                        8, 69, 148]./255;
                end

                figure (1)

                plot(MAct.(['sim_' log_table.Model_Hash{idx_sim_table}]).time,...
                    MAct.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus)*100, ...
                    'LineWidth', line_width,...
                    'Color', sd_colours(i_sd,:), ...
                    'LineStyle',line_style)

                hold on;


            end
        end
        figure(1)
        plot(MAct.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, MAct.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, i_mus)*100, 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')


        xlabel('Movement Duration (%)','FontWeight', 'bold', 'FontSize', axis_label_size);
        ylabel('Muscle Activation (%)' , 'FontWeight', 'bold', 'FontSize', axis_label_size);
        ylim([0, 60]);
        title(modes{i_mode})
        hold on;
        figure(1)
        set(gca,'FontWeight', 'bold')
        %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        saveas(figure(1), [analysis_folder task_name '/Activation_DELT' num2str(i_mus) '_' num2str(modes{i_mode}) '.tiff'],'tiff')
        close all

    end

end



% % % %  %% Analysis muscle moment
% % % % for i_coord = 1:numel(coords_to_plot)
% % % %
% % % %     coord = coords_to_plot{i_coord};
% % % %
% % % %     for i_sim = 1: num_sims
% % % %
% % % %         % Import ForceReporter data
% % % %         muscle_moment = importdata([analysis_folder ['sim_' log_table.Model_Hash{i_sim}] '/Moco_MuscleAnalysis_Moment_' coord '.sto']);
% % % %
% % % %         % Get time
% % % %         MM.(['sim_' log_table.Model_Hash{i_sim}]).time = muscle_moment.data(:,1);
% % % %
% % % %         for i_mus = 1:numel(muscleF_to_plot)
% % % %             % Get handle of coordinates to plot
% % % %             MM.(['sim_' log_table.Model_Hash{i_sim}]).pos(:, i_mus) = find(contains(muscle_moment.colheaders, muscleF_to_plot{i_mus}));
% % % %             MM.(['sim_' log_table.Model_Hash{i_sim}]).label{1,i_mus} = muscle_moment.colheaders{MM.(['sim_' log_table.Model_Hash{i_sim}]).pos(i_mus)};
% % % %
% % % %             MM.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_mus) = muscle_moment.data(:,MM.(['sim_' log_table.Model_Hash{i_sim}]).pos(i_mus));
% % % %
% % % %
% % % %             % Create q plots
% % % %             figure (i_mus)
% % % %
% % % %
% % % %             plot(shoulder_elv_theta.(['sim_' log_table.Model_Hash{i_sim}]), MM.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_mus), 'LineWidth', line_width, 'Color', jrf_colours(i_mus,:))
% % % %
% % % %             hold on;
% % % %             %         scatter(TF.(['sim_' log_table.Model_Hash{i_sim}]).time(max_point + round(numel(TF.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
% % % %             xlabel('Time (s)');
% % % %             ylabel(['Muscle Moment [' muscleF_to_plot{i_mus} '] (Nm)']);
% % % %             hold on;
% % % %         end
% % % %
% % % %     end
% % % %
% % % %
% % % %     % Plot mean scapula data on top of other trials
% % % %     figure(1)
% % % %     plot(MM.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, MM.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 1), 'LineWidth', line_width, 'Color', 'black')
% % % %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % % %     saveas(figure(1), [analysis_folder task_name '/Moment_DELT1_' coord '.tiff'],'tiff')
% % % %     figure(2)
% % % %     plot(MM.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, MM.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 2), 'LineWidth', line_width, 'Color', 'black')
% % % %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % % %     saveas(figure(2), [analysis_folder task_name '/Moment_DELT2_' coord '.tiff'],'tiff')
% % % %     figure(3)
% % % %     plot(MM.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time, MM.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 3), 'LineWidth', line_width, 'Color', 'black')
% % % %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % % %     saveas(figure(3), [analysis_folder task_name '/Moment_DELT3_' coord '.tiff'],'tiff')
% % % %
% % % %     close all
% % % %
% % % %     %% By mode
% % % %     for i_mode = 1:numel(modes)
% % % %
% % % %         for i_sd = 1:numel(SDs)
% % % %
% % % %             morphology = [modes{i_mode} '_' SDs{i_sd}];
% % % %
% % % %             idx_sim_table =  find(contains(log_table.Scapula_morphology, morphology));
% % % %
% % % %             for i_mus = 1:3
% % % %
% % % %                 figure (i_mus)
% % % %                 subplot(2,3, i_mode)
% % % %                 %                 plot(MM.(['sim_' log_table.Model_Hash{idx_sim_table}]).time,...
% % % %                 %                     MM.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus), ...
% % % %                 %                     'LineWidth', line_width,...
% % % %                 %                     'Color', sd_colours(i_sd,:))
% % % %
% % % %                 plot(shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_sim_table}]),...
% % % %                     MM.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus), ...
% % % %                     'LineWidth', line_width,...
% % % %                     'Color', sd_colours(i_sd,:))
% % % %
% % % %                 hold on;
% % % %                 %             scatter(JRA.(['sim_' log_table.Model_Hash{i_sim}]).time(max_point + round(numel(JRA.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
% % % %                 xlabel(['Glenohumeral ABD (' char(176) ')']);
% % % %                 ylabel(['Muscle Moment [' muscleF_to_plot{i_mus} '] (Nm)']);
% % % %                 ylim(y_lims_muscle_moment{i_coord, i_mus});
% % % %                 title(modes{i_mode})
% % % %                 hold on;
% % % %
% % % %             end
% % % %
% % % %
% % % %         end
% % % %
% % % %         figure(1)
% % % %         plot(shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_ave_morph}]), MM.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 1), 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
% % % %         figure(2)
% % % %         plot(shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_ave_morph}]), MM.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 2), 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
% % % %         figure(3)
% % % %         plot(shoulder_elv_theta.(['sim_' log_table.Model_Hash{idx_ave_morph}]), MM.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 3), 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
% % % %
% % % %     end
% % % %
% % % %
% % % %     figure(1)
% % % %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % % %     saveas(figure(1), [analysis_folder task_name '/Moment_vs_ABD_DELT1_' coord '_modes.tiff'],'tiff')
% % % %     figure(2)
% % % %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % % %     saveas(figure(2), [analysis_folder task_name '/Moment_vs_ABD_DELT2_' coord '_modes.tiff'],'tiff')
% % % %     figure(3)
% % % %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% % % %     saveas(figure(3), [analysis_folder task_name '/Moment_vs_ABD_DELT3_' coord '_modes.tiff'],'tiff')
% % % %
% % % %     close all
% % % % end

%% Analysis muscle moment arms
for i_coord = 1:numel(coords_to_plot)

    coord = coords_to_plot{i_coord};

    for i_sim = 1: num_sims

        % Import ForceReporter data
        muscle_moment_arm = importdata([analysis_folder ['sim_' log_table.Model_Hash{i_sim}] '/Moco_MuscleAnalysis_MomentArm_' coord '.sto']);

        % Get time
        MA.(['sim_' log_table.Model_Hash{i_sim}]).time = muscle_moment_arm.data(:,1);

        for i_mus = 1:numel(muscleF_to_plot)



            % Get handle of coordinates to plot
            MA.(['sim_' log_table.Model_Hash{i_sim}]).pos(:, i_mus) = find(contains(muscle_moment_arm.colheaders, muscleF_to_plot{i_mus}));
            MA.(['sim_' log_table.Model_Hash{i_sim}]).label{1,i_mus} = muscle_moment_arm.colheaders{MA.(['sim_' log_table.Model_Hash{i_sim}]).pos(i_mus)};

            %%%%%%%%%%%%%%% To normalise %%%%%%%%%%%%%%%%%%%%%5
            if flag_NormTime == true
                MA.(['sim_' log_table.Model_Hash{i_sim}]).time_og = 100*(muscle_moment_arm.data(:,1)/max(muscle_moment_arm.data(:,1)));
                MA.(['sim_' log_table.Model_Hash{i_sim}]).time = 100*(muscle_moment_arm.data(:,1)/max(muscle_moment_arm.data(:,1)));

                MA.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_mus) = interp1(MA.(['sim_' log_table.Model_Hash{i_sim}]).time_og,...
                    muscle_moment_arm.data(:,MA.(['sim_' log_table.Model_Hash{i_sim}]).pos(i_mus)),...
                    time_vec);
            else
                MA.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_mus) = muscle_moment_arm.data(:,MA.(['sim_' log_table.Model_Hash{i_sim}]).pos(i_mus));

            end

            MA.(['sim_' log_table.Model_Hash{i_sim}]).time = time_vec;
            % Create q plots
            figure (i_mus)


            plot(MA.(['sim_' log_table.Model_Hash{i_sim}]).time,...
                MA.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_mus)*1000, ...
                'LineWidth', line_width, ...
                'Color', jrf_colours(i_mus,:))

            hold on;
            %         scatter(TF.(['sim_' log_table.Model_Hash{i_sim}]).time(max_point + round(numel(TF.(['sim_' log_table.Model_Hash{i_sim}]).data(:, i_joint))*0.6) - 1) , max_F_final(i_sim, i_joint), 'cyan')
            xlabel('Movement Duration (%)','FontWeight', 'bold', 'FontSize', axis_label_size);
            ylabel('Muscle Moment Arm (mm)','FontWeight', 'bold', 'FontSize', axis_label_size);
            hold on;
        end

    end


    % Plot mean scapula data on top of other trials
    figure(1)
    set(gca,'FontWeight', 'bold')
    plot(MA.(['sim_' log_table.Model_Hash{i_sim}]).time, MA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 1)*1000, 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
    %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(1), [analysis_folder task_name '/MomentArm_DELT1_' coord '.tiff'],'tiff')
    figure(2)
    set(gca,'FontWeight', 'bold')
    plot(MA.(['sim_' log_table.Model_Hash{i_sim}]).time, MA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 2)*1000, 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
    %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(2), [analysis_folder task_name '/MomentArm_DELT2_' coord '.tiff'],'tiff')
    figure(3)
    set(gca,'FontWeight', 'bold')
    plot(MA.(['sim_' log_table.Model_Hash{i_sim}]).time, MA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, 3)*1000, 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')
    %     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(figure(3), [analysis_folder task_name '/MomentArm_DELT3_' coord '.tiff'],'tiff')

    save([analysis_folder task_name '/data_MA'], 'MA')
    close all

    %% By mode
    for i_mus = 1:3
        for i_mode = 1:numel(modes)

            for i_sd = 1:numel(SDs)

                morphology = [modes{i_mode} '_' SDs{i_sd}];

                % If comparing two methods of the same mode (i.e. different trials)
                n_idx_sim_table =  find(contains(log_table.Scapula_morphology, morphology));

                for i_trial = 1:numel(n_idx_sim_table)

                    idx_sim_table = n_idx_sim_table(i_trial);

                    if log_table.Date_Generated(idx_sim_table)>=datetime(2023,5,1) && log_table.Date_Generated(idx_sim_table)<=datetime(2023,5,10)
                        line_style = '-';
                        % Reds
                        sd_colours = [254,153,41;
                            236,112,20;
                            204,76,2;
                            140,45,4]./255;
                    elseif log_table.Date_Generated(idx_sim_table)>=datetime(2023,5,10)
                        line_style = '-';
                        % Pinks
                        sd_colours = [221,52,151;
                            174,1,126;
                            122,1,119;
                            73,0,106]./255;
                    else
                        line_style = '-';
                        % Blues
                        sd_colours = [222, 235, 247;
                            158, 202, 225;
                            66, 146, 198;
                            8, 69, 148]./255;
                    end

                    figure (1)
                    plot(MA.(['sim_' log_table.Model_Hash{i_sim}]).time,...
                        MA.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus)*1000, ...
                        'LineWidth', line_width,...
                        'Color', sd_colours(i_sd,:),...
                        'LineStyle', line_style)

                    hold on;

                end
            end
            figure(1)
            plot(MA.(['sim_' log_table.Model_Hash{i_sim}]).time, MA.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, i_mus)*1000, 'LineWidth', line_width, 'Color', 'black','LineStyle','-.')

            xlabel('Movement Duration (%)','FontWeight', 'bold', 'FontSize', axis_label_size);
            ylabel('Muscle Moment Arm (mm)', 'FontWeight', 'bold', 'FontSize', axis_label_size);
            ylim(y_lims_muscle_moment_arm{i_coord, i_mus});
            title(modes{i_mode})
            hold on;

            figure(1)
            set(gca,'FontWeight', 'bold')
            %             set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
            saveas(figure(1), [analysis_folder task_name '/MomentArm_' coord '_DELT' num2str(i_mus) '_' num2str(modes{i_mode}) '.tiff'],'tiff')
            close all


        end

    end

    if strcmp(coord, 'shoulder_elv')
        MA_shoulder_elv = MA;
    elseif strcmp(coord, 'elv_angle')
        MA_elv = MA;
    end
end

for i_mus = 1:3
    for i_mode = 1:numel(modes)

        for i_sd = 1:numel(SDs)

            morphology = [modes{i_mode} '_' SDs{i_sd}];
            % If comparing two methods of the same mode (i.e. different trials)
            n_idx_sim_table =  find(contains(log_table.Scapula_morphology, morphology));

            for i_trial = 1:numel(n_idx_sim_table)

                idx_sim_table = n_idx_sim_table(i_trial);

                if log_table.Date_Generated(idx_sim_table)>=datetime(2023,5,1) && log_table.Date_Generated(idx_sim_table)<=datetime(2023,5,10)
                    line_style = '-';
                    % Reds
                    sd_colours = [254,153,41;
                        236,112,20;
                        204,76,2;
                        140,45,4]./255;

                else
                    line_style = '-';
                    % Blues
                    sd_colours = [222, 235, 247;
                        158, 202, 225;
                        66, 146, 198;
                        8, 69, 148]./255;
                end


                %             figure (1)
                %             plot(MA_elv.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus)*1000,...
                %                 (MA_shoulder_elv.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus)*1000 - MA_shoulder_elv.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(1, i_mus)*1000), ...
                %                 'LineWidth', line_width,...
                %                 'Color', sd_colours(i_sd,:))
                % %  MA_shoulder_elv.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus)*1000 -...
                %             hold on;

                figure(2)
                %             yyaxis left
                subplot(2,1,1)
                plot(MA_elv.(['sim_' log_table.Model_Hash{idx_sim_table}]).time,...
                    MA_shoulder_elv.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus)*1000, ...
                    'LineWidth', line_width,...
                    'LineStyle', '-',...
                    'Color', sd_colours(i_sd,:), ...
                    'Marker', 'none', ...
                    'LineStyle', line_style)
                hold on

                %             yyaxis right
                subplot(2,1,2)
                plot(MA_elv.(['sim_' log_table.Model_Hash{idx_sim_table}]).time,...
                    MA_elv.(['sim_' log_table.Model_Hash{idx_sim_table}]).data(:, i_mus)*1000, ...
                    'LineWidth', line_width,...
                    'LineStyle', '-',...
                    'Color', sd_colours(i_sd,:), ...
                    'Marker', 'none', ...
                    'LineStyle', line_style)
                hold on
            end
        end
        %% Figure one cross hairs
        %         figure(1)
        %         plot(MA_elv.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, i_mus)*1000,...
        %             MA_shoulder_elv.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, i_mus)*1000 - MA_shoulder_elv.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(1, i_mus)*1000,...
        %             'LineWidth', line_width,...
        %             'Color','black', ...
        %             'LineStyle','-')
        %
        %         xlabel(['Muscle Moment Arm Elevation Plane [' muscleF_to_plot{i_mus} '] (mm)'],'FontWeight', 'bold', 'FontSize', axis_label_size);
        %         ylabel(['Muscle Moment Arm Abduction [' muscleF_to_plot{i_mus} '] (mm)'], 'FontWeight', 'bold', 'FontSize', axis_label_size);
        %         %         ylim(y_lims_muscle_moment_arm{1, i_mus});
        %         ylim([-max(abs(y_lims_muscle_moment_arm{1, i_mus})) max(abs(y_lims_muscle_moment_arm{1, i_mus}))]);
        %
        %         xlim([-max(abs(y_lims_muscle_moment_arm{2, i_mus})) max(abs(y_lims_muscle_moment_arm{2, i_mus}))]);
        %
        %         xline(0);
        %         yline(0);
        %         title(modes{i_mode})
        %         hold on;
        %
        %         figure(1)
        %         %             set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        %         saveas(figure(1), [analysis_folder task_name '/MomentArm_ABDvsELV_DELT' num2str(i_mus) '_' num2str(modes{i_mode}) '.tiff'],'tiff')

        %% Figure two Y-axes
        figure(2)

        %         yyaxis left
        subplot(2,1,1)
        plot(MA_elv.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time,...
            MA_shoulder_elv.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, i_mus)*1000, ...
            'LineWidth', line_width,...
            'Color', 'black', ...
            'LineStyle', '-.', ...
            'Marker', 'none')
        %
        %         yyaxis right

        subplot(2,1,2)
        plot(MA_elv.(['sim_' log_table.Model_Hash{idx_ave_morph}]).time,...
            MA_elv.(['sim_' log_table.Model_Hash{idx_ave_morph}]).data(:, i_mus)*1000, ...
            'LineWidth', line_width,...
            'Color', 'black', ...
            'LineStyle', '-.', ...
            'Marker', 'none')
        %
        xlabel('Movement Duration (%)','FontWeight', 'bold', 'FontSize', axis_label_size);
        xticks([0:20:100])
        %         yyaxis left
        subplot(2,1,1)
        ylabel({'Sh. Elv. (mm)'}, 'Color', 'black', 'FontWeight', 'bold', 'FontSize', axis_label_size);
        %         ylim([-10 10]);
        ylim(y_lims_muscle_moment_arm{1, i_mus});
        set(gca,'YColor',[0 0 0], 'FontWeight', 'bold')
        xticks([0:20:100])

        %         yyaxis right
        subplot(2,1,2)
        ylabel({'Elv. Pl. (mm)'}, 'Color', 'black','FontWeight', 'bold', 'FontSize', axis_label_size);
        %         ylim([-10 10])
        ylim(y_lims_muscle_moment_arm{2, i_mus});
        set(gca,'YColor',[0 0 0], 'FontWeight', 'bold')
        xlabel('Movement Duration (%)','FontWeight', 'bold', 'FontSize', axis_label_size);
        xticks([0:20:100])
        hold on;

        figure(2)
        set(gca,'FontWeight', 'bold')
        sgtitle(modes{i_mode},'FontWeight', 'bold')
        %             set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        saveas(figure(2), [analysis_folder task_name '/MomentArm_ABDvsELV2_DELT' num2str(i_mus) '_' num2str(modes{i_mode}) '.tiff'],'tiff')

        close all
    end
end