% shoulder model tests
close all;
clear;
clc;

% Change the current folder to the folder of this m-file.
if ~isdeployed
    cd(fileparts(which(mfilename)));
end

% Import OpenSim 4.3 libraries
import org.opensim.modeling.*

%% Set up

% If replacing Muscles with Actuators
flag_useTorque = false;

% If removing Rotator Cuff muscles
flag_keepRC =  true;

% If altering GH JCS location
flag_changeGHjoint = true;

% If altering GH JCS location and reletive geometry positions of Scapula and Humerus
flag_changeGeom = false;

% Baseline model file
model_name = 'FullShoulderModel_glC.osim';

print_folder_init_str = '1s_con_tol_em1_higherWeight_check_Lck_initGuess_nMusc_500_setTime';

%% Preamble/condition set-up etc.

if flag_useTorque == true

    coords_to_actuate={'shoulder_elv',...
        'shoulder_rot',...
        'elv_angle',...
        'elbow_flexion',...
        'pro_sup'
        };

    actuator_values={30,... 75
        30,... 75
        30,... 30
        30,... 75
        30     %30
        };

    actuator_controls={[1,-1],... % addCoordinateActuator has it as [max,min]
        [1,-1],...
        [1,-1],...
        [1,-1],...
        [1,-1]
        };

end

% Conditions to change JCS position/rotation
if flag_changeGHjoint == true

    % Only for changing position of JCS. 
    % Note: Keep zeros if you don't want to change anything in joint
    conds_GHJ.tran = [ 0 0 0
        ];

    conds_GHJ.rot = [0 0 0;
        ];
elseif flag_changeGeom == true

    % To change JCS position and relative geometry
    conds_GHJ_geom.tran = [%-0.0204 -0.0110 -0.0250 ...  % (1:3) - Position of JCS as defined by first "Descritised" joint <unrothum>
        %-0.01 0.02 -0.02;                               % (4:6) - Position of <Humerus> geometry (relative <Humurus> Body from JCS
        -0.0204 -0.0110 -0.0250 ...
        -0.01 0.02 -0.03;
        ];

    conds_GHJ_geom.rot = [0 0 0 ...
        0 0 0];

end

%% Control tests
% % % start_time=0;
% % % end_time = 2; % seconds
% % % fs = 200; % Sampling frequency (samples per second)
% % %
% % % dt = 1/fs; % seconds per sample
% % %
% % % time_vector = (start_time:dt:end_time)'; % seconds

% Sine wave
% % % F = 1/4; % Sine wave frequency (hertz)

% Ugly but quick (giggity)
% % % sine_control = 35*sin(2*pi*F*time_vector+3*pi/2)+56;
% hold on;plot(time_vector,sine_control)

% Create matrix to pass to OpenSim
% % % coord_data=zeros(size(sine_control,1),22);

% Pass controls to matrix column
% % % coord_data(:,18)=sine_control;

% Create .sto

% % % q.data(:,1)=time_vector;
% % % q.data(:,2:23)=coord_data;
% % %
% % % q.labels={'time'...
% % %     '/jointset/ground_thorax/thorax_tilt/value'...
% % %     '/jointset/ground_thorax/thorax_list/value'...
% % %     '/jointset/ground_thorax/thorax_rotation/value'...
% % %     '/jointset/ground_thorax/thorax_tx/value'...
% % %     '/jointset/ground_thorax/thorax_ty/value'...
% % %     '/jointset/ground_thorax/thorax_tz/value'...
% % %     '/jointset/sternoclavicular/sternoclavicular_r2/value'...
% % %     '/jointset/sternoclavicular/sternoclavicular_r3/value'...
% % %     '/jointset/unrotscap/unrotscap_r3/value'...
% % %     '/jointset/unrotscap/unrotscap_r2/value'...
% % %     '/jointset/acromioclavicular/acromioclavicular_r2/value'...
% % %     '/jointset/acromioclavicular/acromioclavicular_r3/value'...
% % %     '/jointset/acromioclavicular/acromioclavicular_r1/value'...
% % %     '/jointset/unrothum/unrothum_r1/value'...
% % %     '/jointset/unrothum/unrothum_r3/value'...
% % %     '/jointset/unrothum/unrothum_r2/value'...
% % %     '/jointset/shoulder0/elv_angle/value'...
% % %     '/jointset/shoulder1/shoulder_elv/value'...
% % %     '/jointset/shoulder1/shoulder1_r2/value'...
% % %     '/jointset/shoulder2/shoulder_rot/value'...
% % %     '/jointset/elbow/elbow_flexion/value'...
% % %     '/jointset/radioulnar/pro_sup/value'
% % %     };

% % %%  generateMotFile(q.data,q.labels,'coordinates.mot');
%% Modelling

% Check what variable conditions to count to loop correctly
if flag_changeGHjoint == true
    count_var = conds_GHJ;
elseif flag_changeGeom == true
    count_var = conds_GHJ_geom;
end


for i_cond_tran = 1:size(count_var.tran,1)
    for i_cond_rot = 1:size(count_var.rot,1)
        for i_cond_RC = 1:size(flag_keepRC,1)

            osim_model=Model(['..\..\OpenSim\In\Models\' model_name]);
            init_state=osim_model.initSystem();

            ModelVisualizer.addDirToGeometrySearchPaths('..\..\OpenSim\In\Models\Geometry\');

            % Create .sto

            % % % q.data(:,1)=time_vector;
            % % % q.data(:,2:23)=coord_data;

            %%% Create this into a function
            % % % for i_coord=0:osim_model.getCoordinateSet.getSize()-1
            % % %
            % % %     % Get Coordinate name from CoordinateSet of model
            % % %     coord_name=osim_model.getCoordinateSet.get(i_coord).getName();
            % % %
            % % %     % For each coordinate gate Joint
            % % %     joint_name=osim_model.getCoordinateSet.get(i_coord).getJoint().getName();
            % % %
            % % %     % Create char/string of Coordinate to be passed to create .mot/.sto
            % % %     % (/value) by concatination
            % % %     coord_to_write=['/jointset/' char(joint_name) '/' char(coord_name) '/value'];
            % % %
            % % %     disp(coord_to_write);
            % % %
            % % %     % Pass to cell
            % % %     q.labels{i_coord+2}=coord_to_write;
            % % %
            % % %     % Clean up
            % % %     clear coord_name joint_name coord_to_write
            % % %
            % % % end
            % % %
            % % % clear i_coord
            % % %
            % % % % Add time Label
            % % % q.labels{1}='time';
            % % %
            % % % generateMotFile_40(q.data,q.labels,'shoulder_elv.mot');

            % Lock coordinates
            coord_to_lock={'thorax_tilt' 'thorax_list' 'thorax_rotation'...
                'thorax_tx' 'thorax_ty' 'thorax_tz'};

            for i_coord=0:osim_model.getCoordinateSet.getSize()-1

                if nonzeros(strcmp(char(osim_model.getCoordinateSet.get(i_coord).getName()),coord_to_lock))
                    osim_model.getCoordinateSet.get(i_coord).set_locked(true);
                end

            end

            clear i_coord

            %% Add CoordinateActuators (or Reserve) depending if muscles are used in % simulation (torque_flag == false)
            if flag_useTorque==true

                % Remove Muscles and replace with Actuators
                createTorqueDrivenModel(osim_model,...
                    coords_to_actuate,...
                    actuator_values,...
                    actuator_controls,...
                    '_torque');

            else

                % Add coordinate actuators to remaining joints (that don't have muscles)
                for i_coord = 0:osim_model.getCoordinateSet.getSize()-1
                    if contains(char(osim_model.getCoordinateSet.get(i_coord).getName()),'thorax')
                        % Don't add a coordinate actuator as these coordinates won't move
                    elseif strcmp(char(osim_model.getCoordinateSet.get(i_coord).getName()),'elbow_flexion')
                        % Add an idealised torque actuator.
                        addCoordinateActuator(osim_model,char(osim_model.getCoordinateSet.get(i_coord).getName()),...
                            75,...
                            [inf,-inf],...
                            '_torque');

                    elseif strcmp(char(osim_model.getCoordinateSet.get(i_coord).getName()),'pro_sup')
                        % Add an idealised torque actuator.
                        addCoordinateActuator(osim_model,char(osim_model.getCoordinateSet.get(i_coord).getName()),...
                            30,...
                            [inf,-inf],...
                            '_torque');

                    elseif strcmp(char(osim_model.getCoordinateSet.get(i_coord).getName()),'elv_angle')
                        % Add a reserve torque actuator.
                        addCoordinateActuator(osim_model,char(osim_model.getCoordinateSet.get(i_coord).getName()),...
                            1,...
                            [inf,-inf],...
                            '_reserve');

                    elseif strcmp(char(osim_model.getCoordinateSet.get(i_coord).getName()),'shoulder_elv') || ...
                            strcmp(char(osim_model.getCoordinateSet.get(i_coord).getName()),'shoulder_rot')
                        %Add a reserve torque actuator.
                        addCoordinateActuator(osim_model,char(osim_model.getCoordinateSet.get(i_coord).getName()),...
                            1,...
                            [inf,-inf],...
                            '_reserve');

                    end
                end

            end
            clear i_coord
            %% Remove Rotator Cuff muscles (SUPSP, INFSP, SUBSC, TMIN)
            if flag_keepRC(i_cond_RC) == false

                osim_model.updForceSet.remove(osim_model.getMuscles.get('SUPSP'));
                osim_model.updForceSet.remove(osim_model.getMuscles.get('INFSP'));
                osim_model.updForceSet.remove(osim_model.updMuscles.get('SUBSC'));
                osim_model.updForceSet.remove(osim_model.updMuscles.get('TMIN'));

            end

            %% Re-define the GH joints' JCS
            %%% Loop through "Descritised" shoulder model joints
            %%% This methods alters position/orientation of Parent and
            %%% Child (PhysicalOffsetFrames) of Joint but keeps Child
            %%% Body in same ralative position wrt Parent Body.
            if flag_changeGHjoint == true && flag_changeGeom == false

                % "Glenohumeral Joint" (shoulder0, shoulder1, shoulder2)
                % definitions
                joints_to_alter = {'shoulder0',... % 'elv_angle'
                    'shoulder1',...                % 'shoulder_elv' and 'shoulder_rot_2'
                    'shoulder2'};                  % 'shoulder_rot'

                for i_joint = 1:numel(joints_to_alter)

                    %%% JCS tanslation in Parent and Child Body %%%

                    % Get baseline translation in Parent
                    baseline_tran_p = osim_model.getJointSet.get(joints_to_alter{i_joint}).get_frames(0).get_translation();

                    % Update JCS translation in Parent
                    osim_model.updJointSet.get(joints_to_alter{i_joint}).get_frames(0).set_translation(Vec3(...
                        baseline_tran_p.get(0) + conds_GHJ.tran(i_cond_tran, 1),...
                        baseline_tran_p.get(1) + conds_GHJ.tran(i_cond_tran, 2),...
                        baseline_tran_p.get(2) + conds_GHJ.tran(i_cond_tran, 3)));

                    % Get baseline translation in Child
                    baseline_tran_c = osim_model.getJointSet.get(joints_to_alter{i_joint}).get_frames(1).get_translation();

                    % Update JCS translation in Child
                    % Note: Doesn't have to be negated about Z-axis
                    osim_model.updJointSet.get(joints_to_alter{i_joint}).get_frames(1).set_translation(Vec3(...
                        baseline_tran_c.get(0) + conds_GHJ.tran(i_cond_tran, 1),...
                        baseline_tran_c.get(1) + conds_GHJ.tran(i_cond_tran, 2),...
                        baseline_tran_c.get(2) + conds_GHJ.tran(i_cond_tran, 3))); % Note: Negate? Not needed due to SpatialTransform (0 0 -1)

                    %%% JCS orientation/rotation in Parent and Child Body %%%

                    % Get baseline rotation in Parent
                    baseline_rot_p = osim_model.getJointSet.get(joints_to_alter{i_joint}).get_frames(0).get_orientation();

                    % Update JCS rotation in Parent
                    osim_model.updJointSet.get(joints_to_alter{i_joint}).get_frames(0).set_orientation(Vec3(...
                        baseline_rot_p.get(0) + conds_GHJ.rot(i_cond_rot, 1),...
                        baseline_rot_p.get(1) + conds_GHJ.rot(i_cond_rot, 2),...
                        baseline_rot_p.get(2) + conds_GHJ.rot(i_cond_rot, 3)));

                    % Get baseline rotation in Child
                    baseline_rot_c = osim_model.getJointSet.get(joints_to_alter{i_joint}).get_frames(1).get_orientation();

                    % Update JCS rotation in Child
                    osim_model.updJointSet.get(joints_to_alter{i_joint}).get_frames(1).set_orientation(Vec3(...
                        baseline_rot_c.get(0) + conds_GHJ.rot(i_cond_rot, 1),...
                        baseline_rot_c.get(1) + conds_GHJ.rot(i_cond_rot, 2),...
                        baseline_rot_c.get(2) + conds_GHJ.rot(i_cond_rot, 3)));


                end

            elseif flag_changeGHjoint == false && flag_changeGeom == true

                % "Glenohumeral Joint" (shoulder0, shoulder1, shoulder2)
                % definitions
                joints_to_alter = {'unrothum',...
                    'shoulder0',...                % 'elv_angle'
                    'shoulder1',...                % 'shoulder_elv' and 'shoulder_rot_2'
                    'shoulder2'};                  % 'shoulder_rot'

                for i_joint = 1:numel(joints_to_alter)

                    %%% JCS tanslation in Parent and Child Body %%%

                    % Get baseline translation in Parent
                    baseline_tran_p = osim_model.getJointSet.get(joints_to_alter{i_joint}).get_frames(0).get_translation();

                    % Update JCS translation in Parent
                    % Check if it is first joint <unrothum> because that
                    % joint has to define initial position of JCS (for the
                    % following discritised joints [shoulder0/1/2]) and
                    % from there the position of the <Humerus> relative to
                    % the joints' JCS can be defined based on the offset
                    % from the body to Child POF.
                    if strcmp(joints_to_alter{i_joint},'unrothum')
                        osim_model.updJointSet.get(joints_to_alter{i_joint}).get_frames(0).set_translation(Vec3(...
                            baseline_tran_p.get(0) + conds_GHJ_geom.tran(i_cond_tran, 1),...
                            baseline_tran_p.get(1) + conds_GHJ_geom.tran(i_cond_tran, 2),...
                            baseline_tran_p.get(2) + conds_GHJ_geom.tran(i_cond_tran, 3)));
                    else
                        osim_model.updJointSet.get(joints_to_alter{i_joint}).get_frames(0).set_translation(Vec3(...
                            baseline_tran_p.get(0) + conds_GHJ_geom.tran(i_cond_tran, 4),...
                            baseline_tran_p.get(1) + conds_GHJ_geom.tran(i_cond_tran, 5),...
                            baseline_tran_p.get(2) + conds_GHJ_geom.tran(i_cond_tran, 6)));
                    end

                    % Get baseline translation in Child
                    baseline_tran_c = osim_model.getJointSet.get(joints_to_alter{i_joint}).get_frames(1).get_translation();

                    % Update JCS translation in Child
                    % Note: Doesn't have to be negated about Z-axis
                    osim_model.updJointSet.get(joints_to_alter{i_joint}).get_frames(1).set_translation(Vec3(...
                        baseline_tran_c.get(0) + conds_GHJ_geom.tran(i_cond_tran, 4),...
                        baseline_tran_c.get(1) + conds_GHJ_geom.tran(i_cond_tran, 5),...
                        baseline_tran_c.get(2) + conds_GHJ_geom.tran(i_cond_tran, 6))); % Note: Negate? Not needed due to SpatialTransform (0 0 -1)

                    %%% JCS orientation/rotation in Parent and Child Body %%%

                    % Get baseline rotation in Parent
                    baseline_rot_p = osim_model.getJointSet.get(joints_to_alter{i_joint}).get_frames(0).get_orientation();

                    % Update JCS rotation in Parent
                    % Similar check as above Needs to be checked (20220601)
                    if strcmp(joints_to_alter{i_joint},'unrothum')
                        osim_model.updJointSet.get(joints_to_alter{i_joint}).get_frames(0).set_orientation(Vec3(...
                            baseline_rot_p.get(0) + conds_GHJ_geom.rot(i_cond_rot, 1),...
                            baseline_rot_p.get(1) + conds_GHJ_geom.rot(i_cond_rot, 2),...
                            baseline_rot_p.get(2) + conds_GHJ_geom.rot(i_cond_rot, 3)));
                    else
                        osim_model.updJointSet.get(joints_to_alter{i_joint}).get_frames(0).set_orientation(Vec3(...
                            baseline_rot_p.get(0) + conds_GHJ_geom.rot(i_cond_rot, 4),...
                            baseline_rot_p.get(1) + conds_GHJ_geom.rot(i_cond_rot, 5),...
                            baseline_rot_p.get(2) + conds_GHJ_geom.rot(i_cond_rot, 6)));
                    end
                    % Get baseline rotation in Child
                    baseline_rot_c = osim_model.getJointSet.get(joints_to_alter{i_joint}).get_frames(1).get_orientation();

                    % Update JCS rotation in Child
                    osim_model.updJointSet.get(joints_to_alter{i_joint}).get_frames(1).set_orientation(Vec3(...
                        baseline_rot_c.get(0) + conds_GHJ_geom.rot(i_cond_rot, 4),...
                        baseline_rot_c.get(1) + conds_GHJ_geom.rot(i_cond_rot, 5),...
                        baseline_rot_c.get(2) + conds_GHJ_geom.rot(i_cond_rot, 6)));
                end

            end

            % Muscles
            DeGrooteFregly2016Muscle.replaceMuscles(osim_model);

            % Skip this if muscles have been replaced
            if flag_useTorque==false

                % Make problems easier to solve by strengthening the Model and widening the active force-length curve.
                for m = 0:osim_model.getMuscles().getSize()-1
                    musc = osim_model.updMuscles().get(m);
                    musc.setMinControl(0);
                    musc.set_ignore_activation_dynamics(false);
                    musc.set_ignore_tendon_compliance(false);
                    % Strengthens Model
                    musc.set_max_isometric_force(1 * musc.get_max_isometric_force());
                    dgf = DeGrooteFregly2016Muscle.safeDownCast(musc);
                    % Widens Muscle F-L curve
                    dgf.set_active_force_width_scale(1.5);
                    dgf.set_tendon_compliance_dynamics_mode('implicit');
                end

            end

            % Once finished with model changes re-assemble
            osim_model.finalizeConnections();

            % Any other operations on Model should probably be done before it is passed to Problem
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

            % Set qdot Initial and Final
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
            solver.set_optim_max_iterations(5000);

            % % %             % Create an initial guess
            % % %             in_guess=solver.createGuess();
            % % %             in_guess.randomizeAdd();
            % % %             solver.setGuess(in_guess);

            % Set guess from previous (~ good) solution. This will need
            % alterations when other conditions are testeed/
            if flag_keepRC(i_cond_RC) == true
                solver.setGuessFile('..\..\OpenSim\In\Moco\initial_guess\initial_guess_LatReach_RC_1.sto');
            else
                solver.setGuessFile('..\..\OpenSim\In\Moco\initial_guess\initial_guess_LatReach_RC_0.sto')
            end

            %% Solve

            tic;
            predicted_solution = study.solve();
            toc;

            %% Post-processing

            % If failed, unseal
            if ~predicted_solution.success()
                predicted_solution.unseal();
            end

            if flag_changeGHjoint == true

                print_folder_name=horzcat(print_folder_init_str,...
                    '_tran_',...
                    num2str(conds_GHJ.tran(i_cond_tran,1)),...
                    num2str(conds_GHJ.tran(i_cond_tran,2)),...
                    num2str(conds_GHJ.tran(i_cond_tran,3)),...
                    '_rot_',...
                    num2str(conds_GHJ.rot(i_cond_rot,1)),...
                    num2str(conds_GHJ.rot(i_cond_rot,2)),...
                    num2str(conds_GHJ.rot(i_cond_rot,3)),...
                    '_RC_',...
                    num2str(flag_keepRC(i_cond_RC)));

            elseif flag_changeGeom == true

                print_folder_name=horzcat(print_folder_init_str,...
                    '_tran_',...
                    num2str(conds_GHJ_geom.tran(i_cond_tran,1)),...
                    num2str(conds_GHJ_geom.tran(i_cond_tran,2)),...
                    num2str(conds_GHJ_geom.tran(i_cond_tran,3)),...
                    '_rot_',...
                    num2str(conds_GHJ_geom.rot(i_cond_rot,1)),...
                    num2str(conds_GHJ_geom.rot(i_cond_rot,2)),...
                    num2str(conds_GHJ_geom.rot(i_cond_rot,3)),...
                    '_geom_',...
                    num2str(conds_GHJ_geom.tran(i_cond_tran,4)),...
                    num2str(conds_GHJ_geom.tran(i_cond_tran,5)),...
                    num2str(conds_GHJ_geom.tran(i_cond_tran,6)),...
                    '_RC_',...
                    num2str(flag_keepRC(i_cond_RC)));

            end

            if contains (print_folder_name,'.')
                print_folder_name=strrep(print_folder_name,'.','p');
            end

            if contains(print_folder_name, '-')
                print_folder_name=strrep(print_folder_name,'-','m');
            end

            if ~exist(['..\..\OpenSim\Out\Moco\' print_folder_name '\'],"dir")
                mkdir(['..\..\OpenSim\Out\Moco\' print_folder_name '\'])
            end

            new_model_file = ['..\..\OpenSim\Out\Moco\' print_folder_name '\FSModel_GHJoint.osim'];
            solution_file = ['..\..\OpenSim\Out\Moco\' print_folder_name '\MocoSol_LatReach.sto'];

            osim_model.print(new_model_file);
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
            analyzeTool.setModelFilename(new_model_file)
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
    end
end


