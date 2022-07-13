% Parametric modifications of RTSA implant geometry configurations (humeral cup and
% glenosphere) and shoulder anatomy (humerus and scapula) to approximate
% "Virtual Surgery" in RTSA. Main output of this program is to adjust the
% relative position of the humerus with respect to the scapula after the
% RTSA implant components with selected configurations have been positioned
% on the respective anatomies. The configurations include congruent cup and
% hemisphere radii (R_cup = R_hemisphere), position (antero/posterior,
% supero/inferial and base offset) and orientation (antero/postero,
% supero/infero version) of the components. When geometric modification
% have been completed on each anatomy and the implant component positions
% are defined the parametic approximation of the components are exported as
% .stl and then the two anatomies and implant componets are registered in
% the global coordinate system. The newly created shoulder joint CoR in the
% scapula (parent body) and humerus (child body) are exported and redefined
% in OpenSim.
%
% Critical functions:
% axang2rotm, vrrotvec, rotate, fmincon
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% TO ADD  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Reaming depth on glenoid (translate glenoid plane along -Z norm of
%       glenoid plane
% 2) Interactivly chose surgical configuration parameters and updating
%       visualiser and data to pass onto next step?
%       (https://www.mathworks.com/help/control/ug/build-app-with-interactive-plot-updates.html)
% 3) Dynamic reading of .stl files (ALMOST - 20220630)
% 4) Where the hemisphere and cup .stl are saved to 
% 5) Link final GHJ outputs to OpenSim to create new model (COMPLETE
%   20220629) - now to run simulation
% 
%
% Pavlos Silvestros, PhD - University of Victoria, CAN, June 2022

tic
close all;
clear;
clc;

% % % %% Set-up Parallel Computing
% % % 
% % % % Number of Workers
% % % n_workers     = 2;
% % % n_threads      = 36/n_workers;
% % % 
% % % % Specify maximum number of computational threads (?)
% % % maxNumCompThreads(n_threads);
% % % 
% % % % Create parallel pool
% % % pool = parpool('threads', n_workers);
% % % 

%% Define Parameters for hemisphere/cup gemetry and offsets
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hemisphere radius %%%%%%%%%%%%%%%%%%%%%%%%%%%%
diameter = 0.036;

R = diameter/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Glenosphere offsets %%%%%%%%%%%%%%%%%%%%%%%%%%%

hemi_gle_offsets = struct();

% Rotation offsets in degrees

% Anteroversion: +ive; Retroversion: -ive
hemi_gle_offsets.y_ant_retro_version    = 0;
% Inferior inclination: - ive; Superior inclination: +ive
hemi_gle_offsets.x_sup_inf_incl         = -10;

% Translation offsets in meters (m)
hemi_gle_offsets.x_ant_post   = 0;          % X-normal
hemi_gle_offsets.y_prox_dist  = -0.006;     % Y-normal
hemi_gle_offsets.z_base_off   = 0.006;      % Z-normal

%%%%%%%%%%%%%%%%%%%%%%%%%%% Humeral cup offsets %%%%%%%%%%%%%%%%%%%%%%%%%%%

hemi_cup_offsets = struct();

% Rotation offsets in degrees

% Anteroversion: +ive; Retroversion: -ive
hemi_cup_offsets.z_ant_retro_version   = 0;
% Inferior inclination: - ive; Superior inclination: +ive
hemi_cup_offsets.x_sup_inf_incl        = 12.5;

% Translation offsets in meters (m)
hemi_cup_offsets.x_ant_post   = 0;      % X-normal
hemi_cup_offsets.y_base_off   = 0.012;  % Y-normal
hemi_cup_offsets.z_prox_dist  = 0;      % Z-normal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Flags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If should plot intermediate plots for checking
flag_checkPlots = false;

% If replacing Muscles with Actuators
flag_useTorque = false;

% If removing Rotator Cuff muscles
flag_keepRC =  false;

% Replace muscle models Millard2012Equilibrium with DeGrootFregly
flag_ReplaceMuscles = true;

% Run Moco after model is defined?
flag_runSim = true;

% Create a random 11-char hash to reference model file X00yyy111zz (~30e12)
rhash = [char(randi([65 90],1,1))...    
    char(randi([48 57],1,2))...         
    char(randi([97 122],1,3))...       
    char(randi([48 57],1,3))...         
    char(randi([97 122],1,2))];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Call functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create internal functions here (one for humerus and one for scapula) that
% plot and do all the positioning then only return necessary values and
% data.

% Define parametric implant on .stl anatomy & extract parameters in global
scapula = glenoidGeom(R, hemi_gle_offsets, rhash);

% Define parametric implant on .stl anatomy & extract parameters in global
humerus = humerusGeom(R, hemi_cup_offsets, rhash);

% Read in defined implant parameters and .stl and calculate GHJ centre
[GHJ_in_parent, GHJ_in_child] = jointCalculationGH(scapula,humerus);

% Define OpenSim model with new GHJ parameters from 'Virtual Surgery'
model_file = adjustOpenSimModelGHJ(GHJ_in_parent,...
    GHJ_in_child,...
    hemi_gle_offsets,...
    hemi_cup_offsets,...
    R,...
    rhash,...
    flag_useTorque,...
    flag_keepRC,...
    flag_ReplaceMuscles);

% Run OpenSim moco for predictive simulation
if flag_runSim == true
    runRTSAsims(model_file, rhash, flag_keepRC)
end

toc

%% Calculate Joint frame positions with OpenSim definitions
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%% For the humerus
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%% For the scapula
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


