function addTaskBounds(taskName,taskBounds,mocoProblem,modelObject)

% @author: Aaron Fox
% Centre for Sport Research, Deakin University
% aaron.f@deakin.edu.au
% 
% Convenience function for adding bounds to states in a Moco Problem
% relevant to specific tasks.
%
% Input:    taskName - string of the task name being simulated
%           taskBounds - loaded matrix containing task bounds based on Vidt et al. data
%           mocoProblem - MocoProblem object to add goals to
%           modelObject - Opensim model object that is being used in problem
% MODIFIED: Pavlos Silvestros University of Victoria June 2022
    import org.opensim.modeling.*

    %Check for values
    if nargin < 3
        %Throw error
        error('At least 3 inputs (a task name, the task bounds matrix and Moco Problem) are required');
    end

    %% Add bounds.
    %  Note that this function only incorporates the concetric aspect of
    %  the movement.
    
    %Get the index row for the current tasks values in the task
    %bounds structure
    if contains(taskName,'AxillaTouch')
        taskInd = find(strcmp(taskBounds.taskNames,'AxillaTouch'));
    elseif contains(taskName,'ForwardReach')
        taskInd = find(strcmp(taskBounds.taskNames,'ForwardReach'));
    elseif contains(taskName,'HairTouch')
        taskInd = find(strcmp(taskBounds.taskNames,'HairTouch'));
    elseif contains(taskName,'RearTouch')
        taskInd = find(strcmp(taskBounds.taskNames,'RearTouch'));
    elseif contains(taskName,'UpwardReach90')
        taskInd = find(strcmp(taskBounds.taskNames,'UpwardReach90'));
    elseif contains(taskName,'UpwardReach')
        taskInd = find(strcmp(taskBounds.taskNames,'UpwardReach105'));
    elseif contains(taskName,'LateralReach')
        taskInd = find(strcmp(taskBounds.taskNames,'LateralReach'));
    end
    
    %Set relevant task bounds

    %Shoulder elevation
    charValue = ['/jointset/',char(modelObject.getCoordinateSet().get('shoulder_elv').getJoint().getName()),'/',char(modelObject.getCoordinateSet().get('shoulder_elv').getName()),'/value'];
    
    mocoProblem.setStateInfo(charValue,...
        [deg2rad(taskBounds.shoulder_elv.min(taskInd)),...      Range_min
        deg2rad(taskBounds.shoulder_elv.max(taskInd))],...      Range_max
        deg2rad(0),...                                          q_init (can also be a range like below [q_init_min q_init_max])
        [deg2rad(taskBounds.shoulder_elv.con_LB(taskInd)),...   q_fin_min
        deg2rad(taskBounds.shoulder_elv.con_UB(taskInd))]);    %q_fin_max

    %Shoulder rotation
    charValue = ['/jointset/',char(modelObject.getCoordinateSet().get('shoulder_rot').getJoint().getName()),'/',char(modelObject.getCoordinateSet().get('shoulder_rot').getName()),'/value'];
    mocoProblem.setStateInfo(charValue,...
        [deg2rad(taskBounds.shoulder_rot.min(taskInd)),...
        deg2rad(taskBounds.shoulder_rot.max(taskInd))],...
        deg2rad(0),...
        [deg2rad(taskBounds.shoulder_rot.con_LB(taskInd)),...
        deg2rad(taskBounds.shoulder_rot.con_UB(taskInd))]);

    %Elevation plane
    charValue = ['/jointset/',char(modelObject.getCoordinateSet().get('elv_angle').getJoint().getName()),'/',char(modelObject.getCoordinateSet().get('elv_angle').getName()),'/value'];
    mocoProblem.setStateInfo(charValue,...
        [deg2rad(taskBounds.elv_angle.min(taskInd)),...
        deg2rad(taskBounds.elv_angle.max(taskInd))],...
        deg2rad(0),...
        []);
        %[deg2rad(taskBounds.elv_angle.con_LB(taskInd)),...
        %deg2rad(taskBounds.elv_angle.con_UB(taskInd))]);

    %Elbow flexion
    charValue = ['/jointset/',char(modelObject.getCoordinateSet().get('elbow_flexion').getJoint().getName()),'/',char(modelObject.getCoordinateSet().get('elbow_flexion').getName()),'/value'];
    mn = modelObject.getCoordinateSet().get('elbow_flexion').getRangeMin();
    if contains(taskName,'Reach')
        %Limit elbow flexion to 20 degrees
        mx = deg2rad(45);
    else
        %Leave as max attainable elbow flexion
        mx = modelObject.getCoordinateSet().get('elbow_flexion').getRangeMax();
    end
    mocoProblem.setStateInfo(charValue,[mn,mx],0,[])

% % 
% %     mocoProblem.setStateInfo(charValue,...
% %         [deg2rad(0), deg2rad(10)],...
% %         0,...
% %         [deg2rad(0), deg2rad(10)])

    %Forearm
    charValue = ['/jointset/',char(modelObject.getCoordinateSet().get('pro_sup').getJoint().getName()),'/',char(modelObject.getCoordinateSet().get('pro_sup').getName()),'/value'];
    if contains(taskName,'Reach')
        %Limit forearm motion to mostly pronation
        mn = deg2rad(-20);
    else
        %Leave as max attainable supination
        mn = modelObject.getCoordinateSet().get('pro_sup').getRangeMin();
    end
    mx = modelObject.getCoordinateSet().get('pro_sup').getRangeMax();
    mocoProblem.setStateInfo(charValue,[mn,mx],0,[])

%     %Set velocity bounds and for all model coordinates start/end at rest
%     mocoProblem.setStateInfoPattern('/jointset/.*/speed', [-50,50],0,0);
%     
%     %Set muscle activation bounds and for activation to start at zero
%     %(or as close to zero as it should)
%     mocoProblem.setStateInfoPattern('/forceset/.*/activation', [0.01,1],0.01,[]);
    
   
    
    
end