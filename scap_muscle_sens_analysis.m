% scap_muscle_sens_analysis
% Run sensitivity analysis for scapula muscle location definitions
clear
close all
clc

%% Set-up

% Shooulder joint (deg) angle to compute MA for
joint_angle = 90;


% Import OpenSim 4.3 libraries
import org.opensim.modeling.*

% Baseline model file
model_name = 'FSModel_GHJoint_Z87jms113pz.osim';

osim_model=Model(['..\..\OpenSim\In\Models\RTSA_Adjusted\' model_name]);
init_state=osim_model.initSystem();


% Get coordinate handle
shoulder_elv = osim_model.getCoordinateSet().get('shoulder_elv');
elv_angle    = osim_model.getCoordinateSet().get('elv_angle');

%% Update muscle locationtions and (initial) via-point locations
% SSM scapula muscle locations

muscle_names = {'DELT2',...
    'DELT3', ...
    'TRP2',...
    'TRP3', ...
    'TRP4',...
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

        muscle_GP = muscle.getGeometryPath();

        point_name = char(point.getName());

        point_concrete_class = char(point.getConcreteClassName());

        if strcmp(point_concrete_class, 'PathPoint') || strcmp(point_concrete_class, 'ConditionalPathPoint')
            % Downcast AbstractPathPoint to PathPoint to set new xyz
            % location
            point_downCast = PathPoint.safeDownCast(point);

            if strcmp(point_name(end-2:end), 'via') && strcmp('DELT2', char(muscle.getName))
                continue
                %                 % SSM scapula muscle locations - pertubed locations
                %                 muscle_locs = importdata(['..\..\SSM\Scapulas\stl_aligned\m1_0_m2_0_m3_0_' char(muscle.getName) '_pert_data.txt'], ' ');
                %
                %                 % Via point offset from the scap muscle attachement -
                %                 % calculated manually offline
                %                 delt2_via_off = [-0.001, -0.0016, 0.0169];
                %                 location_Vec3 = Vec3(muscle_locs(idx_muscle_data,1) + delt2_via_off(1),...
                %                     muscle_locs(idx_muscle_data,2) + delt2_via_off(2),...
                %                     muscle_locs(idx_muscle_data,3) + delt2_via_off(3));
                %                 point_downCast.set_location(location_Vec3);

            elseif strcmp(point_name(end-2:end), 'via') && strcmp('DELT3', char(muscle.getName))
                continue
                %                 % SSM scapula muscle locations - pertubed locations
                %                 muscle_locs = importdata(['..\..\SSM\Scapulas\stl_aligned\m1_0_m2_0_m3_0_' char(muscle.getName) '_pert_data.txt'], ' ');
                %
                %                 % Via point offset from the scap muscle attachement -
                %                 % calculated manually offline
                %                 delt3_via_off = [-0.0131, -0.0222, 0.0444];
                %                 location_Vec3 = Vec3(muscle_locs(idx_muscle_data,1) + delt3_via_off(1),...
                %                     muscle_locs(idx_muscle_data,2) + delt3_via_off(2),...
                %                     muscle_locs(idx_muscle_data,3) + delt3_via_off(3));
                %                 point_downCast.set_location(location_Vec3);

            elseif strcmp(point_name(end-2:end), 'via') && strcmp('TMAJ', char(muscle.getName))
                continue
                %                 % Via point offset from the scap muscle attachement -
                %                 % calculated manually offline
                %                 tmaj_via_off = [0.0594, 0.0015, 0.0687];
                %                 location_Vec3 = Vec3(muscle_locs(idx_muscle_data,1) + tmaj_via_off(1),...
                %                     muscle_locs(idx_muscle_data,2) + tmaj_via_off(2),...
                %                     muscle_locs(idx_muscle_data,3) + tmaj_via_off(3));
                %                 point_downCast.set_location(location_Vec3);

            else

                % SSM scapula muscle locations - pertubed locations
                muscle_locs = importdata(['..\..\SSM\Scapulas\stl_aligned\m1_0_m2_0_m3_0_' char(muscle.getName) '_pert_data.txt'], ' ');

                for i_pert = 1:size(muscle_locs,1)
                    location_Vec3 = Vec3(muscle_locs(i_pert,1), muscle_locs(i_pert,2), muscle_locs(i_pert,3));
                    point_downCast.set_location(location_Vec3);

                    osim_model.finalizeConnections();
                    osim_model.initSystem();
                    
                    % Set shoulder_elv joint angle
                    osim_model.updCoordinateSet().get('shoulder_elv').setValue(init_state, 30);
                    osim_model.updCoordinateSet().get('elv_angle').setValue(init_state, deg2rad(joint_angle));
                    osim_model.realizePosition(init_state);

                    muscle_MA_perts.(char(muscle.getName()))(i_pert, 1) = muscle_GP.computeMomentArm(init_state, elv_angle);
                end

            end
        else
            error('ERROR: Check path points')
        end
        clear point body_of_point point_downCast location_Vec3
    end
    
    %% Plot moment arms after perubations



    clear muscle muscle_PathPointSet point body_of_point point_downCast location_Vec3
end



for i_muscle = 1:numel(muscle_names)
    scatter(ones(size(muscle_MA_perts.(muscle_names{i_muscle}),1)) + i_muscle -1, muscle_MA_perts.(muscle_names{i_muscle})*1000)
    hold on
end

title(['shoulder_elv MA sensitivity to scapula muscle insertion pertubation @ ' num2str(joint_angle) ' deg'], 'Interpreter', 'none')
xticks(1:numel(muscle_names))
xticklabels(muscle_names)
xlabel("Muscles originating/inserting to the scapula")
ylabel("Muscle moment arm (mm)")

saveas(gcf,['..\..\OpenSim\Out\Verification_Validation\Moment_arms\MA_elv_angle_sensitivity_scap_muscle_position_deg_' num2str(joint_angle) '.fig'], 'fig');
saveas(gcf,['..\..\OpenSim\Out\Verification_Validation\Moment_arms\MA_elv_angle_sensitivity_scap_muscle_position_deg_' num2str(joint_angle) '.tif'], 'tiff');