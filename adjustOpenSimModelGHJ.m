function new_model_file = adjustOpenSimModelGHJ(GHJ_in_parent, GHJ_in_child, hemi_gle_offsets, hemi_cup_offsets, R, rhash, model_SSM, flag_useTorque, flag_keepRC, flag_ReplaceMuscles)
% adjustOpenSimModelGHJ Create new OpenSim model with the newly defined GHJ
% location and with the  generated parametric "implant geometries" 
% (glenosphere and humeral cup). Also export log of RTSA models indexed by 
% 11-charecter alphanumeric hash
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% TO ADD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -
%
% Pavlos Silvestros, PhD - University of Victoria, CAN, June 2022
%% Set-up

% Import OpenSim 4.3 libraries
import org.opensim.modeling.*

% Baseline model file
model_name = 'FullShoulderModel_glC.osim'; %'FullShoulderModelglC_understand_jointDefs_implant.osim';

% To change JCS position and relative geometry
conds_GHJ_geom.tran = [
    GHJ_in_parent ...  % (1, 1:3) - Position of JCS as defined by first "Descritised" joint <unrothum>
    GHJ_in_child       % (1, 4:6) - Position of <Humerus> geometry (relative <Humurus> Body from JCS
    ];

conds_GHJ_geom.rot = [
    0 0 0 ...
    0 0 0];

%% Adjust model

osim_model=Model(['..\..\OpenSim\In\Models\' model_name]);
init_state=osim_model.initSystem();

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
if flag_keepRC == false

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


% "Glenohumeral Joint" (shoulder0, shoulder1, shoulder2)
% definitions
joints_to_alter = {'unrothum',...
    'shoulder0',...                % 'elv_angle'
    'shoulder1',...                % 'shoulder_elv' and 'shoulder_rot_2'
    'shoulder2'};                  % 'shoulder_rot'

for i_joint = 1:numel(joints_to_alter)

    %%% JCS tanslation in Parent and Child Body %%%

    % Get baseline translation in Parent
    %     baseline_tran_p = osim_model.getJointSet.get(joints_to_alter{i_joint}).get_frames(0).get_translation();
    baseline_tran_p = Vec3(0); % [0 0 0];

    % Update JCS translation in Parent
    % Check if it is first joint <unrothum> because that
    % joint has to define initial position of JCS (for the
    % following discritised joints [shoulder0/1/2]) and
    % from there the position of the <Humerus> relative to
    % the joints' JCS can be defined based on the offset
    % from the body to Child POF.
    if strcmp(joints_to_alter{i_joint},'unrothum')
        osim_model.updJointSet.get(joints_to_alter{i_joint}).get_frames(0).set_translation(Vec3(...
            baseline_tran_p.get(0) + conds_GHJ_geom.tran(1, 1),...
            baseline_tran_p.get(1) + conds_GHJ_geom.tran(1, 2),...
            baseline_tran_p.get(2) + conds_GHJ_geom.tran(1, 3)));
    else
        osim_model.updJointSet.get(joints_to_alter{i_joint}).get_frames(0).set_translation(Vec3(...
            baseline_tran_p.get(0) + conds_GHJ_geom.tran(1, 4),...
            baseline_tran_p.get(1) + conds_GHJ_geom.tran(1, 5),...
            baseline_tran_p.get(2) + conds_GHJ_geom.tran(1, 6)));
    end

    % Get baseline translation in Child
    %     baseline_tran_c = osim_model.getJointSet.get(joints_to_alter{i_joint}).get_frames(1).get_translation();
    baseline_tran_c = Vec3(0); %[0 0 0];

    % Update JCS translation in Child
    % Note: Doesn't have to be negated about Z-axis
    osim_model.updJointSet.get(joints_to_alter{i_joint}).get_frames(1).set_translation(Vec3(...
        baseline_tran_c.get(0) + conds_GHJ_geom.tran(1, 4),...
        baseline_tran_c.get(1) + conds_GHJ_geom.tran(1, 5),...
        baseline_tran_c.get(2) + conds_GHJ_geom.tran(1, 6))); % Note: Negate? Not needed due to SpatialTransform (0 0 -1)

    %%% JCS orientation/rotation in Parent and Child Body %%%

    % Get baseline rotation in Parent
    baseline_rot_p = osim_model.getJointSet.get(joints_to_alter{i_joint}).get_frames(0).get_orientation();

    % Update JCS rotation in Parent
    % Similar check as above Needs to be checked (20220601)
    if strcmp(joints_to_alter{i_joint},'unrothum')
        osim_model.updJointSet.get(joints_to_alter{i_joint}).get_frames(0).set_orientation(Vec3(...
            baseline_rot_p.get(0) + conds_GHJ_geom.rot(1, 1),...
            baseline_rot_p.get(1) + conds_GHJ_geom.rot(1, 2),...
            baseline_rot_p.get(2) + conds_GHJ_geom.rot(1, 3)));
    else
        osim_model.updJointSet.get(joints_to_alter{i_joint}).get_frames(0).set_orientation(Vec3(...
            baseline_rot_p.get(0) + conds_GHJ_geom.rot(1, 4),...
            baseline_rot_p.get(1) + conds_GHJ_geom.rot(1, 5),...
            baseline_rot_p.get(2) + conds_GHJ_geom.rot(1, 6)));
    end
    % Get baseline rotation in Child
    baseline_rot_c = osim_model.getJointSet.get(joints_to_alter{i_joint}).get_frames(1).get_orientation();

    % Update JCS rotation in Child
    osim_model.updJointSet.get(joints_to_alter{i_joint}).get_frames(1).set_orientation(Vec3(...
        baseline_rot_c.get(0) + conds_GHJ_geom.rot(1, 4),...
        baseline_rot_c.get(1) + conds_GHJ_geom.rot(1, 5),...
        baseline_rot_c.get(2) + conds_GHJ_geom.rot(1, 6)));
end

% SSM scapula muscle locations
muscle_locs = importdata(['..\..\SSM\Scapulas\stl_aligned\' model_SSM 'muscle_coords.txt'], ' ');
muscle_names = {'DELT2',...
    'DELT3', ...
    'TRAP2',...
    'TRAP3', ...
    'TRAP4',...
    'RMN',...
    'RMJ1',...
    'RMJ2',...
    'CORB',...
    'TMAJ',...
    'PMN',...
    'SRA1',...
    'SRA2',...
    'SRA3',...
    'LVS'};


% Get OpenSim muscle set
muscle_set = osim_model.getMuscles();

% Loop through muscles 
for i_muscle = 0:muscle_set.getSize()-1
    % Get muscle
    muscle = muscle_set.get(i_muscle);
    % Index muscle path point data
    idx_muscle_data = find(strcmp(muscle_names(:), char(muscle.getName)));

    % Skip if muscle isn't of <scapula> 
    if isempty(idx_muscle_data)
        continue
    end
    
    % Get muscle path points
    muscle_PathPointSet = muscle.getGeometryPath().getPathPointSet();
    
    for i_point = 0:muscle_PathPointSet.getSize()-1
        
        point = muscle_PathPointSet.get(i_point);

        body_of_point = point.getBody().getName();
        
        % Skip if point isn' attached to <scapula>
        if ~strcmp(char(body_of_point), 'scapula')
            continue
        end

        % Skip if it is a via point
        point_name = char(point.getName());
        if strcmp(point_name(end-2:end), 'via')
            continue
        end

        point_concrete_class = char(point.getConcreteClassName());
        
        if strcmp(point_concrete_class, 'PathPoint') || strcmp(point_concrete_class, 'ConditionalPathPoint')
            % Downcast AbstractPathPoint to PathPoint to set new xyz
            % location
            point_downCast = PathPoint.safeDownCast(point);
            location_Vec3 = Vec3(muscle_locs(idx_muscle_data,1), muscle_locs(idx_muscle_data,2), muscle_locs(idx_muscle_data,3));
            point_downCast.set_location(location_Vec3);
        else
            error('ERROR: Check path points')
        end

    end
end

% Skip if not replaceing muscle model
if flag_ReplaceMuscles == true

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
end


%% Add/correct meshfile names for the parametric implant geometries


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Scapula %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This will need to add hash or select specific original geometry when we
% analyse different morphologies

% Get the first geometry Scapula
sca_OS_geom = osim_model.getBodySet.get('scapula').get_attached_geometry(0);
sca_OS_mesh = Mesh.safeDownCast(sca_OS_geom);
sca_OS_mesh.set_mesh_file([model_SSM 'downsample.stl'])
% Clone
sca_OS_geom_clone= sca_OS_geom.clone();
% DownCast to <Mesh>
sca_OS_mesh_clone = Mesh.safeDownCast(sca_OS_geom_clone);
% Set new mesh_file_ to clone
sca_OS_mesh_clone.set_mesh_file(string(['gle_' rhash '.stl']));
sca_OS_mesh_clone.setName('scapula_glenosphere');

% Change <Appearance>
sca_OS_mesh_clone.upd_Appearance(0).set_color(Vec3(0,1,0));
sca_OS_mesh_clone.upd_Appearance(0).set_opacity(0.65);

% Attach new Geometry/Mesh to model
osim_model.getBodySet.get('scapula').append_attached_geometry(sca_OS_mesh_clone);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Humerus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the first geometry Humerus
hum_OS_geom = osim_model.getBodySet.get('humerus').get_attached_geometry(0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This will need to add hash or select specific original geometry when we
% analyse different morphologies
% Change name and mesh_fileto resection geometry.
hum_OS_geom.setName('humerus_resected_manifold_closed.stl');
hum_OS_mesh = Mesh.safeDownCast(hum_OS_geom);
hum_OS_mesh.set_mesh_file('humerus_resected_manifold_closed.stl')

% Clone
hum_OS_geom_clone = hum_OS_geom.clone();
% DownCast to <Mesh>
hum_OS_mesh_clone = Mesh.safeDownCast(hum_OS_geom_clone);
% Set new mesh_file_ to clone
hum_OS_mesh_clone.setName('humeral_cup');
hum_OS_mesh_clone.set_mesh_file(string(['cup_' rhash '.stl']));

% Change <Appearance>
hum_OS_mesh_clone.upd_Appearance(0).set_color(Vec3(1,0,0));
hum_OS_mesh_clone.upd_Appearance(0).set_opacity(0.5);

% Attach new Geometry/Mesh to model
osim_model.getBodySet.get('humerus').append_attached_geometry(hum_OS_mesh_clone);

% Once finished with model changes re-assemble
osim_model.finalizeConnections();



%% Print this bad boy

% Add random pause between 0.25-0.50 seconds to print files in parfor
pause(0.250 + rand*0.250)

date_time_now = datestr(datetime);

f_id = fopen('..\..\OpenSim\In\Models\RTSA_Adjusted\RTSA_model_log.txt', 'a+');


fprintf(f_id, '#########################################################\n\n');
fprintf(f_id, '%s produced at %s \r\n', rhash, date_time_now);
fprintf(f_id, ['\nModel RTSA configuration (values in m and degrees): \r\n' ...
    '\nScapula morphology (SSM model):\t %s \r\n' ...
    '\nRTSA joint position: \r\n' ...
    ['\tGHJ in Parent (Scapula):\t' repmat('%g ', 1, numel(GHJ_in_parent)), '\r\n'] ...
    ['\tGHJ in Child (Humerus):\t\t' repmat('%g ', 1, numel(GHJ_in_child)), '\r\n']...
    '\nGlenosphere Configuration: \r\n'...
    '\tAntero/Retro version:\t\t%g \r\n' ...
    '\tSupero/Infero incl/tion:\t%g \r\n' ...
    '\tBaseplate offset:\t\t\t%g \r\n' ...
    '\tAntero/Poserior offset:\t\t%g \r\n' ...
    '\tSupero/Inferior offset:\t\t%g \r\n' ...
    '\tRadius:\t\t\t\t%g \r\n' ...
    '\nHumeral cup Configuration: \r\n'...
    '\tAntero/Retro version:\t\t%g \r\n' ...
    '\tSupero/Infero incl/tion:\t%g \r\n' ...
    '\tBaseplate offset:\t\t\t%g \r\n' ...
    '\tAntero/Poserior offset:\t\t%g \r\n' ...
    '\tSupero/Inferior offset:\t\t%g \r\n' ...
    '\tRadius:\t\t\t\t%g \r\n\r\n'], ...
    string(model_SSM), ...
    GHJ_in_parent, ...
    GHJ_in_child, ...
    hemi_gle_offsets.y_ant_retro_version, ...
    hemi_gle_offsets.x_sup_inf_incl, ...
    hemi_gle_offsets.z_base_off, ...
    hemi_gle_offsets.x_ant_post, ...
    hemi_gle_offsets.y_prox_dist, ...
    R,...
    hemi_cup_offsets.z_ant_retro_version, ...
    hemi_cup_offsets.x_sup_inf_incl, ...
    hemi_cup_offsets.y_base_off, ...
    hemi_cup_offsets.x_ant_post, ...
    hemi_cup_offsets.z_prox_dist, ...
    R);

fclose(f_id);

% Read individual model table and .txt in logs
var_names_table = {'Date_Generated',...
    'Model_Hash',...
    'Scapula_morphology', ...
    'GHJ_in_parent',...
    'GHJ_in_child',...
    'hemi_gle_offsets',...
    'hemi_radius',...
    'hemi_cup_offsets',...
    'cup_radius'};

txt_table = splitvars(table(string(date_time_now),...
    string(rhash),...
    string(model_SSM),...
    GHJ_in_parent,...
    GHJ_in_child, ...
    struct2table(hemi_gle_offsets),...
    R,...
    struct2table(hemi_cup_offsets),...
    R,...
    'VariableNames',var_names_table));

writetable(txt_table, '..\..\OpenSim\In\Models\RTSA_Adjusted\RTSA_model_log_table.txt');

% Write OpenSim model
new_model_file = ['..\..\OpenSim\In\Models\RTSA_Adjusted\FSModel_GHJoint_' rhash '.osim'];
osim_model.print(new_model_file);

disp('OpenSim model with adjusted GHJ definition and geomtry for RTSA has been defined and saved at: ')
disp(['..\..\OpenSim\In\Models\RTSA_Adjusted\FSModel_GHJoint_' rhash '.osim'])

disp('GHJ joint position (x y z) in Parent (Scapula) is: ')
disp(GHJ_in_parent)

disp('GHJ joint position (x y z) in Child (Humerus) is: ')
disp(GHJ_in_child)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here a function that pases the model to the simulation pipeline


