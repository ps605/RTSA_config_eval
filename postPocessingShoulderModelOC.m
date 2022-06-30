function postPocessingShoulderModelOC(solution_file, print_folder_name, jra_filename)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

coords_to_plot = {'elv_angle',...
    'shoulder_elv',...
    'shoulder_rot',...
    'elbow_flexion'};

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
    102 194 164
    44 127 184;
    221 28 119]./255;

jointF_to_plot = {'shoulder0_on_scapphant_in_glenoid_centre_fx'...
    'shoulder0_on_scapphant_in_glenoid_centre_fy',...
    'shoulder0_on_scapphant_in_glenoid_centre_fz'};

%% Coordinates (q)
states = importdata(solution_file);

% Get time
q.time = states.data(:,1);
mus.time = q.time;

%%% Handle q data
for i_coord = 1:numel(coords_to_plot)
    % Get handle of coordinates to plot
    q.pos(:, i_coord) = find(contains(states.colheaders, [coords_to_plot{i_coord} '/value']));
    q.label{1,i_coord} = states.colheaders{q.pos(i_coord)};
    q.data(:, i_coord) = rad2deg(states.data(:,q.pos(i_coord)));

    % Create q plots
    figure (10)
    plot(q.time, q.data(:, i_coord), 'LineWidth',1.5)
    xlabel('Time (s)');
    ylabel('q (deg)');
    hold on;

end

% Add legend
legend(coords_to_plot,...
    'Location','northeastoutside',...
    'Interpreter', 'none',...
    'AutoUpdate', 'off');


% Save plot
saveas(figure(10),['..\..\OpenSim\Out\Moco\' print_folder_name '\coordinates.tif'], 'tif' )
saveas(figure(10),['..\..\OpenSim\Out\Moco\' print_folder_name '\coordinates.fig'], 'fig' )
close;


%% Activations (a)
%%% Handle activation data

% Check if Rotator Cuffs exist
if contains(solution_file, 'RC_0')
    for i_mus = 5:numel(muscles_to_plot)

        % Get handle of mus to plot
        mus.pos(:, i_mus) = find(contains(states.colheaders, [muscles_to_plot{i_mus} '/activation']));
        mus.label{1,i_mus} = states.colheaders{mus.pos(i_mus)};
        mus.data(:, i_mus) = states.data(:,mus.pos(i_mus));

        % Create mus plots
        figure (11)
        plot(mus.time, mus.data(:, i_mus), 'LineWidth', 1.5, 'Color', muscle_colours(i_mus,:))
        xlabel('Time (s)');
        ylabel('Activation [0-1]');
        ylim([0 1]);
        hold on;

    end


    % Add legend
    legend(muscles_to_plot(5:end),...
        'Location','northeastoutside',...
        'Interpreter', 'none',...
        'AutoUpdate', 'off');
else
    for i_mus = 1:numel(muscles_to_plot)

        % Get handle of mus to plot
        mus.pos(:, i_mus) = find(contains(states.colheaders, [muscles_to_plot{i_mus} '/activation']));
        mus.label{1,i_mus} = states.colheaders{mus.pos(i_mus)};
        mus.data(:, i_mus) = states.data(:,mus.pos(i_mus));

        % Create mus plots
        figure (11)
        plot(mus.time, mus.data(:, i_mus), 'LineWidth',1.5, 'Color', muscle_colours(i_mus,:))
        xlabel('Time (s)');
        xlabel('Time (s)');
        ylabel('Activation [0-1]');
        ylim([0 1]);
        hold on;

    end


    % Add legend
    legend(muscles_to_plot,...
        'Location','northeastoutside',...
        'Interpreter', 'none',...
        'AutoUpdate', 'off');
end


saveas(figure(11),['..\..\OpenSim\Out\Moco\' print_folder_name '\activations.tiff'], 'tif' )
saveas(figure(11),['..\..\OpenSim\Out\Moco\' print_folder_name '\activations.fig'], 'fig' )

close;

%% JRA

joint_reaction = importdata(jra_filename);


% Get time
JRA.time = joint_reaction.data(:,1);


%%% Handle q data
for i_joint = 1:numel(jointF_to_plot)
    % Get handle of coordinates to plot
    JRA.pos(:, i_joint) = find(contains(joint_reaction.colheaders, jointF_to_plot{i_joint}));
    JRA.label{1,i_joint} = joint_reaction.colheaders{JRA.pos(i_joint)};

    if contains(JRA.label{1,i_joint}, 'fz')
        JRA.data(:, i_joint) = abs(joint_reaction.data(:,JRA.pos(i_joint)));
    else
        JRA.data(:, i_joint) = joint_reaction.data(:,JRA.pos(i_joint));
    end

    % Create q plots
    figure (12)
    plot(JRA.time, JRA.data(:, i_joint), 'LineWidth', 1.5, 'Color', jrf_colours(i_joint,:))
    xlabel('Time (s)');
    ylabel('JRF (N)');
    hold on;

end

% Plot Resultant Load
JRA.resultant = sqrt(JRA.data(:, 1).^2 + JRA.data(:, 2).^2 + JRA.data(:, 3).^2);

figure (12)
plot(JRA.time, JRA.resultant(:, 1), 'LineWidth', 1.5, 'Color', jrf_colours(4,:))
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


% Save plot
saveas(figure(12),['..\..\OpenSim\Out\Moco\' print_folder_name '\JRA.tif'], 'tif' )
saveas(figure(12),['..\..\OpenSim\Out\Moco\' print_folder_name '\JRA.fig'], 'fig' )
close;




end