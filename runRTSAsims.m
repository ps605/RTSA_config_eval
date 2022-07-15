function runRTSAsims(model_file, rhash, flag_keepRC)
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
problem.setTimeBounds(0,0.550);

% Set coordinate q bounds
% coords_to_bound={'elv_angle',...
%     'shoulder_elv',...
%     'shoulder_rot',...
%     'elbow_flexion',...
%     'pro_sup'};
%
% coord_bound={[-0.0242407, 1.36556]
%     [0, 1.85581]
%     [-0.138833, 0.713998]
%     [0, 1.5708]
%     [-1.5708, 1.5708]};
%
% for i_coord=1:numel(coords_to_bound)
%
%     % Coord name
%     coord_char_val=['/jointset/',...
%     char(osim_model.getCoordinateSet.get(coords_to_bound{i_coord}).getJoint.getName()),...
%     '/',...
%     char(osim_model.getCoordinateSet.get(coords_to_bound{i_coord}).getName()),...
%     '/value'];
%
%     % Coord bound
%    % coord_bound=[osim_model.getCoordinateSet.get(coords_to_bound{i_coord}).getRangeMin(), osim_model.getCoordinateSet.get(coords_to_bound{i_coord}).getRangeMax()];
%
%     problem.setStateInfoPattern(coord_char_val, coord_bound{i_coord}, 0, []);
% end

% This should be changed back to task_bounds - the ones used in Fox for first
% proper run
load('new_task_bounds.mat');
addTaskBounds('LateralReach',new_task_bounds,problem,osim_model);

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

%problem.addGoal(MocoFinalTimeGoal('time',1));

%% 2 - Minimise distance from defined point at final time (From Aaron Fox, Deakin)

% % % %Get the position of the C7 marker
% % % model_state = osim_model.initSystem();
% % % C7 = osim_model.getMarkerSet().get('C7').getLocationInGround(model_state);
% % %
% % % %Prescribe the marker end point using the X and Z coordinates of
% % % %the C7 marker and add the arbitrary distance to the Y position
% % % hairReachPoint = Vec3(C7.get(0),C7.get(1) + 0.25,C7.get(2));
% % %
% % % % Lateral reach point
% % %

% % % %Cleanup
% % % clear C7

% % %             % Add the end point costs with appropriate weights
% % %             endPointCost1 = MocoMarkerFinalGoal('MF_endPoint',5);
% % %             endPointCost1.setPointName('/markerset/MiddleFinger');
% % %
% % %             % Point at at approximatly shoulder hieght.
% % %             lat_reach_point=Vec3(0.05, -0.021, 0.85);
% % %             endPointCost1.setReferenceLocation(lat_reach_point);
% % %
% % %             %Add the end point cost along with an effort cost.
% % %             problem.addGoal(endPointCost1);

% Three Marker end point costs (all in Ground)
point_RS = Vec3(0.0796176, -0.0852458, 0.675867);
point_US = Vec3(0.015706, -0.0795034, 0.687447);
point_wri_out = Vec3(0.0508631, -0.056026, 0.685693);

end_point_cost_1 = MocoMarkerFinalGoal('marker_RS', 60);
end_point_cost_1.setPointName('/markerset/RS');
end_point_cost_1.setReferenceLocation(point_RS);


end_point_cost_2 = MocoMarkerFinalGoal('marker_US', 60);
end_point_cost_2.setPointName('/markerset/US');
end_point_cost_2.setReferenceLocation(point_US);


end_point_cost_3 = MocoMarkerFinalGoal('marker_wrist_out', 60);
end_point_cost_3.setPointName('/markerset/wri_out');
end_point_cost_3.setReferenceLocation(point_wri_out);

% Add the MarkerGoals to the problem
problem.addGoal(end_point_cost_1);
problem.addGoal(end_point_cost_2);
problem.addGoal(end_point_cost_3);

%% 3 - Minimise effort

problem.addGoal(MocoControlGoal('effort',10));

%% Configure the solver.

solver = study.initCasADiSolver();
solver.set_num_mesh_intervals(75);
solver.set_verbosity(2);
solver.set_optim_solver('ipopt');
solver.set_optim_convergence_tolerance(1e-1);
solver.set_optim_constraint_tolerance(1e-3);
solver.set_optim_max_iterations(1);

% % %             % Create an initial guess
% % %             in_guess=solver.createGuess();
% % %             in_guess.randomizeAdd();
% % %             solver.setGuess(in_guess);

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

% If failed, unseal
if ~predicted_solution.success()
    predicted_solution.unseal();
end

print_folder_name = ['sim_' rhash];

if ~exist(['..\..\OpenSim\Out\Moco\' print_folder_name '\'],"dir")
    mkdir(['..\..\OpenSim\Out\Moco\' print_folder_name '\'])
end

solution_file = ['..\..\OpenSim\Out\Moco\' print_folder_name '\MocoSol_LatReach.sto'];

predicted_solution.write(['..\..\OpenSim\Out\Moco\' print_folder_name '\MocoSol_LatReach.sto']);

%% ANALYSIS

states_storage=Storage(solution_file);

% Set ArrayStr
joints = ArrayStr();
joints.set(0, 'shoulder0');

on_bodies = ArrayStr();
on_bodies.set(0, 'parent');

in_frames = ArrayStr();
in_frames.set(0, '/bodyset/scapula/glenoid_centre/');

% Get time
time_array = ArrayDouble;
states_storage.getTimeColumn(time_array);

% Set up AnalyzeTool
analyzeTool=AnalyzeTool('..\..\OpenSim\In\Setup_files\JRA\templateJRA.xml',0);
analyzeTool.setName('Moco');
analyzeTool.setInitialTime(0);
analyzeTool.setFinalTime(time_array.getLast);
%             analyzeTool.setStatesStorage(states_storage);
analyzeTool.setStatesFileName(solution_file);
%             analyzeTool.setModel(osim_model);
analyzeTool.setModelFilename(model_file)
analyzeTool.setResultsDir(['..\..\OpenSim\Out\Moco\' print_folder_name '\']);

JRA = analyzeTool.updAnalysisSet.get(0);
%             JRA.setName('ShoulderJointReaction')
%             JRA.setJointNames(joints);
%             JRA.setOnBody(on_bodies)
%             JRA.setInFrame(in_frames);
JRA.setEndTime(time_array.getLast)

% Print and read back in as bug workaround
analyzeTool.print('runAnalyzeTool.xml');

runTool = AnalyzeTool('runAnalyzeTool.xml');
runTool.run();

%% POST PROCESSING

jra_filename = ['..\..\OpenSim\Out\Moco\' print_folder_name '\Moco_JointReaction_ReactionLoads.sto'];
postPocessingShoulderModelOC(solution_file,print_folder_name, jra_filename);


end