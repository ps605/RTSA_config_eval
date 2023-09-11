function addCoordinateActuator(modelObject,coordinate,optForce,controlLevel,appendStr)

% @author: Aaron Fox
% Centre for Sport Research, Deakin University
% aaron.f@deakin.edu.au
% 
% Convenience function for adding a coordinate actuator to model
%
% Input:    modelObject - Opensim model object to add actuator to
%           coordinate - string of coordinate name for actuator
%           optForce - value for actuators optimal force
%           controlLevel - [x,y] values for max and min control
%           appendStr - string to append to coordinate name in setting actuator name

    import org.opensim.modeling.*

    %Check for values
    if nargin < 5
        appendStr = '_torque';
    end
    
    if nargin < 4
        controlLevel = [inf,-inf];
    end
    
    if nargin < 3
        %Throw error
        error('At least 3 inputs (a model object, coordinate name and optimal force value) are required');
    end

    %% Create actuator
    actu = CoordinateActuator();
    
    %Set name
    actu.setName([coordinate,appendStr]);
    
    %Set coordinate for actuator
    actu.setCoordinate(modelObject.getCoordinateSet().get(coordinate));
    
    %Set optimal force
    actu.setOptimalForce(optForce);
    
    %Set control levels    
    actu.setMaxControl(controlLevel(1));
    actu.setMinControl(controlLevel(2));
    
    %Add actuator to model
    modelObject.addComponent(actu);
    
end