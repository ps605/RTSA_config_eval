function runRTSAsims(model_file, rhash, flag_keepRC, task_name)
%runRTSAsims Summary of this function goes here
%   Detailed explanation goes here


%% Set-up

% Import OpenSim 4.3 libraries
import org.opensim.modeling.*

osim_model = Model(model_file);

%% Simulation

% Moco
study=MocoStudy();

% Initialise and Access optimisation Problem
problem=study.updProblem();

% Pass Model (system dynamics) to Problem
problem.setModel(osim_model);

%% Bounds and constraints

% Set time bounds
problem.setTimeBounds(0,[0.75 1.5]);


% This should be changed back to task_bounds - the ones used in Fox for first
% proper run
load('new_task_bounds.mat');
addTaskBounds(task_name,new_task_bounds,problem,osim_model);

clear i_coord

% Set q_dot Initial and Final
problem.setStateInfoPattern('/jointset/.*/speed', [-50,50],0,0);

% Set muscle activation bounds
problem.setStateInfoPattern('/forceset/.*/activation', [0.01,1],0.01,[]);

% Explicity constraint locked coordinates
coord_to_lock={'thorax_tilt' 'thorax_list' 'thorax_rotation'...
    'thorax_tx' 'thorax_ty' 'thorax_tz'};


for i_coord=1:numel(coord_to_lock)

    problem.setStateInfo(['/jointset/ground_thorax/' coord_to_lock{i_coord} '/value'], [0,0],0,0);
    problem.setStateInfo(['/jointset/ground_thorax/' coord_to_lock{i_coord} '/speed'], [0,0],0,0);

end

clear i_coord
%% Goals
%% 1 - Minimise time taken

problem.addGoal(MocoFinalTimeGoal('time',1));

%% 2 - Minimise distance from defined point at final time (From Aaron Fox, Deakin)

if strcmp(task_name, 'HairTouch')
    % Get the position of the C7 marker
    init_state = osim_model.initSystem();
    C7 = osim_model.getMarkerSet().get('C7').getLocationInGround(init_state);

    %Prescribe the marker end point using the X and Z coordinates of
    %the C7 marker and add the arbitrary distance to the Y position
    point_occiput = Vec3(C7.get(0),C7.get(1) + 0.25,C7.get(2));

    end_point_cost_1 = MocoMarkerFinalGoal('marker_MiddleFinger', 60);
    end_point_cost_1.setPointName('/markerset/MiddleFinger');
    end_point_cost_1.setReferenceLocation(point_occiput);

    % Add the MarkerGoals to the problem
    problem.addGoal(end_point_cost_1);

elseif strcmp(task_name, 'LateralReach')

    init_state = osim_model.initSystem();

    SJC_ground = osim_model.getJointSet().get('shoulder0').get_frames(1).getPositionInGround(init_state);

    %Calculate the distance of the forearm (i.e. between the elbow and wrist
    %joint centre).

    %Get the position of the joint centres. Joint 1 corresponds to ulna offset
    %frame for elbow and joint 0 the radius offset for the radius hand joint
    EJC_ground = osim_model.getJointSet().get('elbow').get_frames(1).getPositionInGround(init_state);
    WJC_ground = osim_model.getJointSet().get('radius_hand_r').get_frames(0).getPositionInGround(init_state);

    %Calculate the distance between the joint centres
    elbow = [EJC_ground.get(0),EJC_ground.get(1),EJC_ground.get(2)];
    wrist = [WJC_ground.get(0),WJC_ground.get(1),WJC_ground.get(2)];
    FA_length = dist_markers(elbow,wrist);

    %Calculate the position 200% forearm length in front of the shoulder. In
    %front is represented by positive X
    point_lateral_SJC = [SJC_ground.get(0), SJC_ground.get(1), SJC_ground.get(2)+(FA_length*2.25)];

    RS_marker = osim_model.getMarkerSet().get('RS').getLocationInGround(init_state);
    US_marker = osim_model.getMarkerSet().get('US').getLocationInGround(init_state);
    RS = [RS_marker.get(0),RS_marker.get(1),RS_marker.get(2)];
    US = [US_marker.get(0),US_marker.get(1),US_marker.get(2)];
    wristWidth = dist_markers(RS,US);

    %Add and subtract half of the wrist distance from the original marker end
    %point along the X-axis to get the proposed end points for the markers. It
    %is positive X in the ground frame for the radius marker and negative X for
    %the ulna marker
    point_US = Vec3(point_lateral_SJC(1)-(wristWidth/2), point_lateral_SJC(2), point_lateral_SJC(3));
    point_RS = Vec3(point_lateral_SJC(1)+(wristWidth/2), point_lateral_SJC(2), point_lateral_SJC(3));

    %Measure the distance from the wrist joint centre to the wri_out marker for
    %prescribing where the hand needs to go.
    wri_out = osim_model.getMarkerSet().get('wri_out').getLocationInGround(init_state);
    wri_out = [wri_out.get(0),wri_out.get(1),wri_out.get(2)];
    wrist = [WJC_ground.get(0),WJC_ground.get(1),WJC_ground.get(2)];
    wristHeight = dist_markers(wri_out,wrist);

    %Add the wirst height amount along the y-axis from the proposed reach point
    %to get the point where the wri_out marker needs to go
    point_wri_out = Vec3(point_lateral_SJC(1),point_lateral_SJC(2)+wristHeight,point_lateral_SJC(3));

    % % Three Marker end point costs (all in Ground)
    % point_RS = Vec3(0.0796176, -0.0852458, 0.675867);
    % point_US = Vec3(0.015706, -0.0795034, 0.687447);
    % point_wri_out = Vec3(0.0508631, -0.056026, 0.685693);

    end_point_cost_1 = MocoMarkerFinalGoal('marker_RS', 60);
    end_point_cost_1.setPointName('/markerset/RS');
    end_point_cost_1.setReferenceLocation(point_RS);
    % Define end_point marker postion for visualisation
    osim_model.getMarkerSet().get('LateralREach_RS').set_location(point_RS)


    end_point_cost_2 = MocoMarkerFinalGoal('marker_US', 60);
    end_point_cost_2.setPointName('/markerset/US');
    end_point_cost_2.setReferenceLocation(point_US);
    % Define end_point marker postion for visualisation
    osim_model.getMarkerSet().get('LateralReach_US').set_location(point_US);

    end_point_cost_3 = MocoMarkerFinalGoal('marker_wrist_out', 60);
    end_point_cost_3.setPointName('/markerset/wri_out');
    end_point_cost_3.setReferenceLocation(point_wri_out);
    % Define end_point marker postion for visualisation
    osim_model.getMarkerSet().get('LateralReach_wri_out').set_location(point_wri_out);

    % Add the MarkerGoals to the problem
    problem.addGoal(end_point_cost_1);
    problem.addGoal(end_point_cost_2);
    problem.addGoal(end_point_cost_3);

elseif strcmp(task_name, 'UpwardReach')

    init_state = osim_model.initSystem();

    SJC_ground = osim_model.getJointSet().get('shoulder0').get_frames(1).getPositionInGround(init_state);

    %Calculate the distance of the forearm (i.e. between the elbow and wrist
    %joint centre).

    %Get the position of the joint centres. Joint 1 corresponds to ulna offset
    %frame for elbow and joint 0 the radius offset for the radius hand joint
    EJC_ground = osim_model.getJointSet().get('elbow').get_frames(1).getPositionInGround(init_state);
    WJC_ground = osim_model.getJointSet().get('radius_hand_r').get_frames(0).getPositionInGround(init_state);

    %Calculate the distance between the joint centres
    elbow = [EJC_ground.get(0),EJC_ground.get(1),EJC_ground.get(2)];
    wrist = [WJC_ground.get(0),WJC_ground.get(1),WJC_ground.get(2)];
    FA_length = dist_markers(elbow,wrist);

    %Calculate the position 200% forearm length in front of the shoulder. In
    %front is represented by positive X
    point_anterior_SJC = [SJC_ground.get(0)+(FA_length*2.0), SJC_ground.get(1), SJC_ground.get(2)];

    RS_marker = osim_model.getMarkerSet().get('RS').getLocationInGround(init_state);
    US_marker = osim_model.getMarkerSet().get('US').getLocationInGround(init_state);
    RS = [RS_marker.get(0),RS_marker.get(1),RS_marker.get(2)];
    US = [US_marker.get(0),US_marker.get(1),US_marker.get(2)];
    wristWidth = dist_markers(RS,US);

    %Add and subtract half of the wrist distance from the original marker end
    %point along the Z-axis to get the proposed end points for the markers. It
    %is positive Z in the ground frame for the ulnar marker and negative Z for
    %the radius marker
    point_US = Vec3(point_anterior_SJC(1), point_anterior_SJC(2), point_anterior_SJC(3)+(wristWidth/2));
    point_RS = Vec3(point_anterior_SJC(1), point_anterior_SJC(2), point_anterior_SJC(3)-(wristWidth/2));

    %Measure the distance from the wrist joint centre to the wri_out marker for
    %prescribing where the hand needs to go.
    wri_out = osim_model.getMarkerSet().get('wri_out').getLocationInGround(init_state);
    wri_out = [wri_out.get(0),wri_out.get(1),wri_out.get(2)];
    wrist = [WJC_ground.get(0),WJC_ground.get(1),WJC_ground.get(2)];
    wristHeight = dist_markers(wri_out,wrist);

    %Add the wirst height amount along the y-axis from the proposed reach point
    %to get the point where the wri_out marker needs to go
    point_wri_out = Vec3(point_anterior_SJC(1),point_anterior_SJC(2)+wristHeight,point_anterior_SJC(3));

    % % Three Marker end point costs (all in Ground)
    % point_RS = Vec3(0.0796176, -0.0852458, 0.675867);
    % point_US = Vec3(0.015706, -0.0795034, 0.687447);
    % point_wri_out = Vec3(0.0508631, -0.056026, 0.685693);

    end_point_cost_1 = MocoMarkerFinalGoal('marker_RS', 60);
    end_point_cost_1.setPointName('/markerset/RS');
    end_point_cost_1.setReferenceLocation(point_RS);
    % Define end_point marker postion for visualisation
    osim_model.getMarkerSet().get('LateralREach_RS').set_location(point_RS)

    end_point_cost_2 = MocoMarkerFinalGoal('marker_US', 60);
    end_point_cost_2.setPointName('/markerset/US');
    end_point_cost_2.setReferenceLocation(point_US);
    % Define end_point marker postion for visualisation
    osim_model.getMarkerSet().get('LateralReach_US').set_location(point_US);

    end_point_cost_3 = MocoMarkerFinalGoal('marker_wrist_out', 60);
    end_point_cost_3.setPointName('/markerset/wri_out');
    end_point_cost_3.setReferenceLocation(point_wri_out);
    % Define end_point marker postion for visualisation
    osim_model.getMarkerSet().get('LateralReach_wri_out').set_location(point_wri_out);

    % Add the MarkerGoals to the problem
    problem.addGoal(end_point_cost_1);
    problem.addGoal(end_point_cost_2);
    problem.addGoal(end_point_cost_3);

end

%% 3 - Minimise effort

problem.addGoal(MocoControlGoal('effort',10));

%% Configure the solver.

solver = study.initCasADiSolver();
solver.set_num_mesh_intervals(75);
solver.set_verbosity(2);
solver.set_optim_solver('ipopt');
solver.set_optim_convergence_tolerance(1e-1);
solver.set_optim_constraint_tolerance(1e-3);
solver.set_optim_max_iterations(10);

% % % % Create an initial guess
% % % in_guess=solver.createGuess();
% % % in_guess.randomizeAdd();
% % % solver.setGuess(in_guess);

% Set guess from previous (~ good) solution. This will need
% alterations when other conditions are testeed/
if flag_keepRC == true
    solver.setGuessFile('..\..\OpenSim\In\Moco\initial_guess\initial_guess_LatReach_RC_1.sto');
else
    solver.setGuessFile('..\..\OpenSim\In\Moco\initial_guess\initial_guess_LatReach_RC_0.sto')
end

%% Solve

tic;
predicted_solution = study.solve();
toc;

% Add random pause between 0.25-1 seconds to print files in parfor
pause(0.250 + rand*0.075)
%% Post-processing
% Finilise connections for redefined marker targets and print model 
osim_model.finalizeConnections();
osim_model.print(model_file);

% If failed, unseal
if ~predicted_solution.success()
    predicted_solution.unseal();
end

print_folder_name = ['sim_' rhash];

if ~exist(['..\..\OpenSim\Out\Moco\' print_folder_name '\'],"dir")
    mkdir(['..\..\OpenSim\Out\Moco\' print_folder_name '\'])
end

solution_file = ['..\..\OpenSim\Out\Moco\' print_folder_name '\MocoSol_' task_name '.sto'];

predicted_solution.write(['..\..\OpenSim\Out\Moco\' print_folder_name '\MocoSol_' task_name '.sto']);

%% ANALYSIS

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
analyzeTool.setResultsDir(['..\..\OpenSim\Out\Moco\' print_folder_name '\']);

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

%% POST PROCESSING

jra_filename = ['..\..\OpenSim\Out\Moco\' print_folder_name '\Moco_JointReaction_ReactionLoads.sto'];
postPocessingShoulderModelOC(solution_file, print_folder_name, jra_filename);


end