% postProcessCompare

clear ;
close all;
clc;

% Change the current folder to the folder of this m-file.
if ~isdeployed
    cd(fileparts(which(mfilename)));
end


orthoload.mean = 607.60;
orthoload.std = 142.28;

% What batch of sims to look for in ...\Out\Moco\
lookup_sim_string = 'nMusc_3Point_LatRch_MI_100*';

clean_labels = {'Native Shoulder',...
    'Native Shoulder - No RC',...
    'No RC - Medial = 10 mm',...
    'No RC - Medial = 20 mm',...
    'No RC - Inferior = 10 mm',...
    'No RC - Inferior = 10 mm and Medial 20 mm'
    };

% List the folder with data to compare
list_sims_to_compare = dir(['..\..\OpenSim\Out\Moco\' lookup_sim_string]);

% Sellect which trials to plo
openvar('list_sims_to_compare');
keyboard

% JRF variables to plot
jointF_to_plot = {'shoulder0_on_scapphant_in_glenoid_centre_fx'...
    'shoulder0_on_scapphant_in_glenoid_centre_fy',...
    'shoulder0_on_scapphant_in_glenoid_centre_fz'};

jrf_colours = [240 59 32;
    102 194 164
    44 127 184;
    221 28 119]./255;

% Activations to plot
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



for i_sim = 1 : numel(list_sims_to_compare)

    solution_file = ['..\..\OpenSim\Out\Moco\' list_sims_to_compare(i_sim).name '\MocoSol_LatReach.sto'];
    jra_filename = ['..\..\OpenSim\Out\Moco\' list_sims_to_compare(i_sim).name '\' list_sims_to_compare(i_sim).name '_JointReaction_ReactionLoads.sto'];

    %% JRA

    % Import data
    joint_reaction = importdata(jra_filename);

    % Get sim name
    JRA.sim_name{i_sim,1} = list_sims_to_compare(i_sim).name;

    % Identify underscore locations for indexing conditions
    underscore_locs=strfind(JRA.sim_name(i_sim),'_');
    
    JRA.sim_label{i_sim,1} = strrep(JRA.sim_name{i_sim,1}(underscore_locs{1}(end-5)+1:end),'_',' ');

    % Get time
    JRA.time(:,1,i_sim) = joint_reaction.data(:,1);


    %%% Handle q data
    for i_joint = 1:numel(jointF_to_plot)
        % Get handle of coordinates to plot
        JRA.pos(:, i_joint,i_sim) = find(contains(joint_reaction.colheaders, jointF_to_plot{i_joint}));
        JRA.label{1,i_joint,i_sim} = joint_reaction.colheaders{JRA.pos(:, i_joint,i_sim)};

        if contains(JRA.label{1,i_joint,i_sim}, 'fz')
            JRA.F_data(:, i_joint,i_sim) = abs(joint_reaction.data(:,JRA.pos(:, i_joint,i_sim)));
        else
            JRA.F_data(:, i_joint,i_sim) = joint_reaction.data(:,JRA.pos(:, i_joint,i_sim));
        end

        % Create q plots

        figure (12)

        if contains(list_sims_to_compare(i_sim).name, 'RC_1')
            plot(JRA.time(:,1,i_sim), JRA.F_data(:, i_joint,i_sim), '-', 'LineWidth', 0.5, 'Color', jrf_colours(i_joint,:))
        elseif contains(list_sims_to_compare(i_sim).name, 'RC_0')
            plot(JRA.time(:,1,i_sim), JRA.F_data(:, i_joint,i_sim), '--', 'LineWidth', 0.5, 'Color', jrf_colours(i_joint,:))
        end

        xlabel('Time (s)');
        ylabel('JRF (N)');
        hold on;

    end

    % Instantiate max Field
    JRA.F_max(i_sim,1:4) = zeros;

    % Get max values for x-y-z components
    % In X
    [max_Fx_v, max_Fx_p] = max(abs(JRA.F_data(:, 1,i_sim)));
    JRA.F_max(i_sim,1) = JRA.F_data(max_Fx_p, 1,i_sim);

    % In Y
    [max_Fy_v, max_Fy_p] = max(abs(JRA.F_data(:, 2,i_sim)));
    JRA.F_max(i_sim,2) = JRA.F_data(max_Fy_p, 2,i_sim);

    % In Z
    [max_Fz_v, max_Fz_p] = max(abs(JRA.F_data(:, 3,i_sim)));
    JRA.F_max(i_sim,3) = JRA.F_data(max_Fz_p, 3,i_sim);

%     JRA.F_max(i_sim,1:3) = max(JRA.F_data(:, :,i_sim),[],1);

    % Calculate Resultant Load
    JRA.F_res(:,1,i_sim) = sqrt(JRA.F_data(:, 1,i_sim).^2 + JRA.F_data(:, 2,i_sim).^2 + JRA.F_data(:, 3,i_sim).^2);

    % Get max for F_res
    JRA.F_max(i_sim,4) = max(JRA.F_res(:,1,i_sim),[],1);

    % Plot Resultant Vector
    figure (12)

    if contains(list_sims_to_compare(i_sim).name, 'RC_1')
        plot(JRA.time(:,1,i_sim), JRA.F_res(:,1,i_sim), '-', 'LineWidth', 2, 'Color', jrf_colours(4,:))
    elseif contains(list_sims_to_compare(i_sim).name, 'RC_0')
        plot(JRA.time(:,1,i_sim), JRA.F_res(:,1,i_sim), '--', 'LineWidth', 2, 'Color', jrf_colours(4,:))
    end
    xlabel('Time (s)');
    ylabel('JRF (N)');

    % Add legend
    legend({'Glenoid Anterior (+) - Posterior (-) Shear',...
        'Glenoid Superior (+) - Inferior (-) Shear',...
        'Glenoid Compression',...
        'Resultant Glenoid Joint Load'},...
        'Location','southoutside',...
        'Interpreter', 'none',...
        'AutoUpdate', 'off');
end

% Plot Orthoload data ma
yline(orthoload.mean, '-', 'LineWidth', 2);
yline(orthoload.mean + orthoload.std, '-', 'LineWidth', 1);
yline(orthoload.mean - orthoload.std, '-', 'LineWidth', 1);

% Add legend
legend({'Glenoid Anterior (+) / Posterior (-) Shear',...
    'Glenoid Superior (+) / Inferior (-) Shear',...
    'Glenoid Compression',...
    'Resultant Glenoid Joint Load'},...
    'Location','southoutside',...
    'Interpreter', 'none',...
    'AutoUpdate', 'off');


% Plot max values
h = figure(13);

plot_labels = categorical(clean_labels); % JRA.sim_label
plot_labels = reordercats(plot_labels,clean_labels); % JRA.sim_label

b = bar(plot_labels, JRA.F_max);

% Plot Orthoload data ma
yline(orthoload.mean, '-', 'LineWidth', 2.5);
yline(orthoload.mean + orthoload.std, '-', 'LineWidth', 1.5);
yline(orthoload.mean - orthoload.std, '-', 'LineWidth', 1.5);

ylabel('JRF Max (N)');
xlabel('Shoulder Model Conditions')
title('Lateral Reach Task - Joint Reaction Forces')
for i_c = 1:numel(b)
    b(i_c).FaceColor = jrf_colours(i_c,:);
end



ylim([-500 2000])

legend({'Glenoid Anterior (+) / Posterior (-) Shear',...
    'Glenoid Superior (+) / Inferior (-) Shear',...
    'Glenoid Compression',...
    'Resultant Glenoid Joint Load'},...
    'Location','northeast',...
    'Interpreter', 'none',...
    'AutoUpdate', 'off');

% Make some hard changes to plot
set(figure(13),'Position',[0 365  690   630])
h.Children(2).FontSize = 11;
h.Children(2).FontWeight = 'bold';
h.Children(2).TickDir = 'none';


% Save plot
mkdir(['..\..\OpenSim\Out\Moco\Analysis\' lookup_sim_string(1:end-1) '\']);

saveas(figure(13),['..\..\OpenSim\Out\Moco\Analysis\' lookup_sim_string(1:end-1) '\JRA_comp.tif'], 'tif' )
saveas(figure(13),['..\..\OpenSim\Out\Moco\Analysis\' lookup_sim_string(1:end-1) '\JRA_comp.fig'], 'fig' )
close;
