close all;
clear;
clc;

% Import OpenSim 4.3 libraries
import org.opensim.modeling.*

model_name='FullShoulderModel.osim';

osim_model=Model(['..\..\OpenSim\In\Models\' model_name]);
init_state=osim_model.initSystem();

% safeDownCast Body::<scapula> to PhysicalFrame
PF_scapula=PhysicalFrame.safeDownCast(osim_model.getBodySet.get('scapula'));

% Create Frame
POF_glenoid_centre = PhysicalOffsetFrame();
POF_glenoid_centre.setName('glenoid_centre');
POF_glenoid_centre.setParentFrame(PF_scapula);

% Centre of Glenoid fossa @ scapular in FullShoulderModel.osim
POF_glenoid_centre.set_translation(Vec3(-0.0300675, -0.0391136, -0.015)); 

% Rotate for x-axis as parallel as posible to glenoid surf
POF_glenoid_centre.set_orientation(Vec3(0, 0.55, 0)); % @ scapular in FullShoulderModel.osim

% Pass back to the Model
osim_model.updBodySet.get('scapula').addComponent(POF_glenoid_centre);

% Finalise Connections
osim_model.finalizeConnections;

osim_model.print(['..\..\OpenSim\In\Models\' model_name(1:end-5) 'glC.osim']);