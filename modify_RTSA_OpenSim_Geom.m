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
%
% Pavlos Silvestros, PhD - University of Victoria, CAN, June 2022

close all;
clear;
clc;

% Time the run
time_i = datetime;
%% Set-up

%%%%%%%%%%%%%%%%% Create parameter combinations for loops %%%%%%%%%%%%%%%%%
design_param.diameter       = {0.036};
design_param.hemi_base_off  = {0};

scapula_morphologies = {'m1_0_m2_0_m3_0_'};

% Create permutation matrix
param_matrix= allcomb( ...
    design_param.diameter,...
    design_param.hemi_base_off, ...
    scapula_morphologies...
    );

% Split matrix
p1 = param_matrix(:,1);
p2 = param_matrix(:,2);
p3 = param_matrix(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Flags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If use parallel with 2 workers (18 threads/worker) to batch job
flag_useParallel = false;

% If should plot intermediate plots for checking
flag_checkPlots = false;

% If replacing Muscles with Actuators
flag_useTorque = false;

% If removing Rotator Cuff muscles
flag_keepRC =  false;

% Replace muscle models Millard2012Equilibrium with DeGrootFregly
flag_ReplaceMuscles = true;

% Run Moco after model is defined?
flag_runSim = false;

% Optimise DELT1, DELT2 and DELT3 via points
flag_viaPointOpt = true;
flag_DELT1 = false;
flag_DELT2 = true;
flag_DELT3 = false;
%% Pass setup parameters and prepare models/simulations
if flag_useParallel == true

    %%%%%%%%%%%%%%%%%%%%%%%%% Parallel Computing %%%%%%%%%%%%%%%%%%%%%%%%%%

    % Number of Workers
    n_workers     = 2;
    n_threads      = 36/n_workers;

    % Specify maximum number of computational threads (?)
    % maxNumCompThreads(n_threads);

    % Create parallel pool
    pool = parpool('2Workers');

    parfor i_param = 1:size(param_matrix,1)
        %% Define Parameters for hemisphere/cup gemetry and offsets
       
        %%%%%%%%%%%%%%%%%%%%%%% Hemisphere radius %%%%%%%%%%%%%%%%%%%%%%%%%
        diameter = p1{i_param};

        R = diameter/2;

        %%%%%%%%%%%%%%%%%%%%%%% Glenosphere offsets %%%%%%%%%%%%%%%%%%%%%%%

        hemi_gle_offsets = struct();

        % Rotation offsets in degrees

        % Anteroversion: +ive; Retroversion: -ive
        hemi_gle_offsets.y_ant_retro_version    = 0;
        % Inferior inclination: - ive; Superior inclination: +ive
        hemi_gle_offsets.x_sup_inf_incl         = -10;

        % Translation offsets in meters (m)
        hemi_gle_offsets.x_ant_post   = 0;          % X-normal
        hemi_gle_offsets.y_prox_dist  = -0.006;     % Y-normal
        hemi_gle_offsets.z_base_off   = p2{i_param};      % Z-normal

        %%%%%%%%%%%%%%%%%%%%%%% Humeral cup offsets %%%%%%%%%%%%%%%%%%%%%%%

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

        %%%%%%%%%%%%%%%%%%%%%%%% Model morphology %%%%%%%%%%%%%%%%%%%%%%%%%
        model_SSM = p3{i_param};

        % Create a random 11-char hash to reference model file X00yyy111zz (~30e12)
        % Add random pause between 0.25-0.50 seconds to print files in parfor
        pause(0.250 + rand*0.250)
        rng('shuffle');
        rhash = [char(randi([65 90],1,1))...
            char(randi([48 57],1,2))...
            char(randi([97 122],1,3))...
            char(randi([48 57],1,3))...
            char(randi([97 122],1,2))];

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%% Call functions %%%%%%%%%%%%%%%%%%%%%%%%%
        % Create internal functions here (one for humerus and one for scapula) that
        % plot and do all the positioning then only return necessary values and
        % data.

        % Define parametric implant on .stl anatomy & extract parameters in global
        scapula = glenoidGeom(R, hemi_gle_offsets, model_SSM, rhash);

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

        close all

        % Run OpenSim moco for predictive simulation
        if flag_runSim == true
            runRTSAsims(model_file, rhash, flag_keepRC)
        end

    end

elseif flag_useParallel == false
    for i_param = 1:size(param_matrix,1)
        %% Define Parameters for hemisphere/cup gemetry and offsets
   
        %%%%%%%%%%%%%%%%%%%%%%%% Hemisphere radius %%%%%%%%%%%%%%%%%%%%%%%%
        diameter = p1{i_param};

        R = diameter/2;

        %%%%%%%%%%%%%%%%%%%%%%% Glenosphere offsets %%%%%%%%%%%%%%%%%%%%%%%

        hemi_gle_offsets = struct();

        % Rotation offsets in degrees

        % Anteroversion: +ive; Retroversion: -ive
        hemi_gle_offsets.y_ant_retro_version    = 0;
        % Inferior inclination: - ive; Superior inclination: +ive
        hemi_gle_offsets.x_sup_inf_incl         = -10;

        % Translation offsets in meters (m)
        hemi_gle_offsets.x_ant_post   = 0;          % X-normal
        hemi_gle_offsets.y_prox_dist  = -0.006;     % Y-normal
        hemi_gle_offsets.z_base_off   = p2{i_param};      % Z-normal

        %%%%%%%%%%%%%%%%%%%%%%% Humeral cup offsets %%%%%%%%%%%%%%%%%%%%%%%

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
        
        %%%%%%%%%%%%%%%%%%%%%%%% Model morphology %%%%%%%%%%%%%%%%%%%%%%%%%
        model_SSM = p3{i_param};

        % Create a random 11-char hash to reference model file X00yyy111zz (~30e12)
        % Add random pause between 0.25-0.50 seconds to print files in parfor
        pause(0.250 + rand*0.250)
        rng('shuffle');
        rhash = [char(randi([65 90],1,1))...
            char(randi([48 57],1,2))...
            char(randi([97 122],1,3))...
            char(randi([48 57],1,3))...
            char(randi([97 122],1,2))];

        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%% Call functions %%%%%%%%%%%%%%%%%%%%%%%%%
        % Create internal functions here (one for humerus and one for scapula) that
        % plot and do all the positioning then only return necessary values and
        % data.

        % Define parametric implant on .stl anatomy & extract parameters in global
        scapula = glenoidGeom(R, hemi_gle_offsets, model_SSM, rhash);

        % Define parametric implant on .stl anatomy & extract parameters in global
        humerus = humerusGeom(R, hemi_cup_offsets, rhash);

        % Read in defined implant parameters and .stl and calculate GHJ centre
        [GHJ_in_parent, GHJ_in_child] = jointCalculationGH(scapula,humerus);

        close 10 20 2

        % Define OpenSim model with new GHJ parameters from 'Virtual Surgery'
        model_file = adjustOpenSimModelGHJ(GHJ_in_parent,...
            GHJ_in_child,...
            hemi_gle_offsets,...
            hemi_cup_offsets,...
            R,...
            rhash,...
            model_SSM,...
            flag_useTorque,...
            flag_keepRC,...
            flag_ReplaceMuscles);
        
        if flag_viaPointOpt  == true
            optimDeltViaPoint(model_file, flag_DELT1, flag_DELT2, flag_DELT3)
        end
        
        % Run OpenSim moco for predictive simulation
        if flag_runSim == true
            runRTSAsims(model_file, rhash, flag_keepRC)
        end

    end
end

% Show entire time of simulation batch
time_f = datetime;
run_time = time_f - time_i;

disp('#######################################################');
disp('Overall simulation batch took....');
disp(' ');
disp(run_time);
disp('#######################################################');
