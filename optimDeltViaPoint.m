function optimDeltViaPoint(model_file, flag_DELT1, flag_DELT2, flag_DELT3)
%optimDeltViaPoint Defines DELT1, DELT2 and DELT3 via poit locations
%   Detailed explanation goes here

%% Set-up
import org.opensim.modeling.*

osim_model = Model(model_file);
init_state = osim_model.initSystem();

% Ackland et al (2010) RTSA MA data

data_RTSA.angles        = [2.5, 30, 60, 90, 120];
data_RTSA.DELT1_mean    = [15.6, 25.2,32.5, 35.8, 33.3]*0.001;
data_RTSA.DELT1_sd      = [2.3, 2.5, 1.8, 3.4, 2.7]*0.001;
data_RTSA.DELT2_mean    = [30.2, 33.9, 42.2, 46.2, 39.8]*0.001;
data_RTSA.DELT2_sd      = [6.6, 5.5, 5.5, 4.2, 5.8]*0.001;
data_RTSA.DELT3_mean    = [1.3, 3.5, 7.3, 11.4, 14.1]*0.001;
data_RTSA.DELT3_sd      = [1.4, 1.0, 1.7, 2.5, 3.6]*0.001;

%% Handle model

% Get muscle handles
delt1 = osim_model.getMuscles.get('DELT1');
delt2 = osim_model.getMuscles.get('DELT2');
delt3 = osim_model.getMuscles.get('DELT3');

% Get coordinate handle
shoulder_elv = osim_model.getCoordinateSet().get('shoulder_elv');

% Get PathPointSet(s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELT1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get GeometryPath to calculate MomentArm later
delt1_GP = delt1.getGeometryPath();
delt1_PPS = delt1_GP.getPathPointSet();
delt1_PPS_size = delt1_PPS.getSize();

% Loop through PathPointSet to get via point
for i_point = 0:delt1_PPS_size-1

    point = delt1_PPS.get(i_point);
    point_name = char(point.getName());

    if strcmp(point_name(end-2:end), 'via')

        % Get via point and downcast to get/set location
        delt1_via_downCast = PathPoint.safeDownCast(point);

        % Get via point location as vector for handling- not Vec3
        delt1_via_loc = [delt1_via_downCast.get_location().get(0),...
            delt1_via_downCast.get_location().get(1),...
            delt1_via_downCast.get_location().get(2)];

    else
        continue
    end
end

% DELT1 moment arm before any optimisation
for i_angle = 1:length(data_RTSA.angles)
    osim_model.updCoordinateSet().get('shoulder_elv').setValue(init_state, deg2rad(data_RTSA.angles(i_angle)));
    osim_model.realizePosition(init_state);

    model_MA_init.DELT1(i_angle) = delt1_GP.computeMomentArm(init_state,shoulder_elv);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELT2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get GeometryPath to calculate MomentArm later
delt2_GP = delt2.getGeometryPath();
delt2_PPS = delt2_GP.getPathPointSet();
delt2_PPS_size = delt2_PPS.getSize();

% Some check for it being via point HERE
% Loop through PathPointSet to get via point
for i_point = 0:delt2_PPS_size-1

    point = delt2_PPS.get(i_point);
    point_name = char(point.getName());

    if strcmp(point_name(end-2:end), 'via')

        % Get via point and downcast to get/set location
        delt2_via_downCast = PathPoint.safeDownCast(point);

        % Get via point location as vector for handling- not Vec3
        delt2_via_loc = [delt2_via_downCast.get_location().get(0),...
            delt2_via_downCast.get_location().get(1),...
            delt2_via_downCast.get_location().get(2)];

    else
        continue
    end
end

% DELT2 moment arm before any optimisation
for i_angle = 1:length(data_RTSA.angles)
    osim_model.updCoordinateSet().get('shoulder_elv').setValue(init_state, deg2rad(data_RTSA.angles(i_angle)));
    osim_model.realizePosition(init_state);

    model_MA_init.DELT2(i_angle) = delt2_GP.computeMomentArm(init_state,shoulder_elv);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELT3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get GeometryPath to calculate MomentArm later
delt3_GP = delt3.getGeometryPath();
delt3_PPS = delt3_GP.getPathPointSet();
delt3_PPS_size = delt3_PPS.getSize();

for i_point = 0:delt3_PPS_size-1

    point = delt3_PPS.get(i_point);
    point_name = char(point.getName());

    if strcmp(point_name(end-2:end), 'via')

        % Get via point and downcast to get/set location
        delt3_via_downCast = PathPoint.safeDownCast(point);

        % Get via point location as vector for handling- not Vec3
        delt3_via_loc = [delt3_via_downCast.get_location().get(0),...
            delt3_via_downCast.get_location().get(1),...
            delt3_via_downCast.get_location().get(2)];

    else
        continue
    end
end

% DELT3 moment arm before any optimisation
for i_angle = 1:length(data_RTSA.angles)
    osim_model.updCoordinateSet().get('shoulder_elv').setValue(init_state, deg2rad(data_RTSA.angles(i_angle)));
    osim_model.realizePosition(init_state);

    model_MA_init.DELT3(i_angle) = delt3_GP.computeMomentArm(init_state,shoulder_elv);
end


%% Optimisation - MA calculation after via point and joint modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELT1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_DELT1 == true
    % Search radius around init location
    radius = 0.025;
    p_sim_0 = delt1_via_loc;

    ub = p_sim_0 + radius;% delt1_via_loc + radius;%[0.05, 0.05, 0.05];
    lb = p_sim_0 - radius; %delt1_via_loc - radius;%[-0.05, -0.05, -0.05];
    lb(3) = p_sim_0(3);

    figure(101);
    hold on
    title('via point location in x-y-z space (m)')
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    scatter3(p_sim_0(1), p_sim_0(2), p_sim_0(3), 'o', 'filled','magenta')


    % Inequality constraint to keep new via point location with radius of x of
    % original via point
    %%% fCon = @(p_sim)sphere_func_con(p_sim, p_sim_0, radius);
    % Cost function to minimise moment arm differances between simulated and
    % calculated conditions
    fObj = @(p_sim)J_momentArmDist(p_sim, data_RTSA, osim_model, 'DELT1', delt1_via_downCast);

    % Set-up options
    options = optimoptions('ga', 'Display', 'iter', 'PlotFcn',{@gaplotbestf, @gaplotmaxconstr});
    options.FunctionTolerance       = 1e-6;
    options.ConstraintTolerance     = 1e-6;
    options.MaxGenerations          = 30;
    options.PopulationSize          = 40;
    options.EliteCount              = ceil(0.05*options.PopulationSize);
    options.FitnessLimit            = 1e-5;
    options.InitialPopulationMatrix = []; % [0.0016, -0.0163, 0.0148]; %[0.0176, -0.0090, 0.0210];%[0.0271, 0.0048, 0.0189]; %[0.0272, 0.0048, 0.0189]; %[-0.0258, 0.0189, 0.0198];
    options.CreationFcn             = [];
    options.MaxStallGenerations     = 15;
    % options.MaxFunctionEvaluations  = 1000;
    % options.StepTolerance           = 1e-20;
    % options.FiniteDifferenceStepSize = 0.001;

    % Create problem structure to back to optimiser
    problem.fitnessfcn   =   fObj;
    problem.nvars       =   3;
    problem.options     =   options;
    problem.x0          =   p_sim_0;
    problem.Aineq       =   [];
    problem.bineq       =   [];
    problem.Aeq         =   [];
    problem.beq         =   [];
    problem.lb          =   lb;
    problem.ub          =   ub;
    problem.nonlcon     =   [];
    problem.solver      =   'ga';

    % Run fmincon
    [p_sim_opt,J_opt,~,~] = ga(problem);

    % Optimised point for DELT1
    opt_via(1,1:3) = p_sim_opt;
    clear('p_sim_opt')

    plotViaOptResults('DELT1',...
        p_sim_0,...
        opt_via(1,1:3),...
        delt1_via_downCast,...
        delt1_GP,...
        shoulder_elv,...
        osim_model,...
        data_RTSA,...
        model_MA_init.DELT1,...
        J_opt,...
        radius);




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELT2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif flag_DELT2 == true
    % Search radius around init location
    radius = 0.025;
    p_sim_0 = delt2_via_loc;

    ub = p_sim_0 + radius;% delt1_via_loc + radius;%[0.05, 0.05, 0.05];
    lb = p_sim_0 - radius; %delt1_via_loc - radius;%[-0.05, -0.05, -0.05];
    lb(3) = p_sim_0(3);

    figure(101);
    hold on
    title('via point location in x-y-z space (m)')
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    scatter3(p_sim_0(1), p_sim_0(2), p_sim_0(3), 'o', 'filled','magenta')


    % Inequality constraint to keep new via point location with radius of x of
    % original via point
    %%% fCon = @(p_sim)sphere_func_con(p_sim, p_sim_0, radius);
    % Cost function to minimise moment arm differances between simulated and
    % calculated conditions
    fObj = @(p_sim)J_momentArmDist(p_sim, data_RTSA, osim_model, 'DELT2', delt2_via_downCast);

    % Set-up options
    options = optimoptions('ga', 'Display', 'iter', 'PlotFcn',{@gaplotbestf, @gaplotmaxconstr});
    options.FunctionTolerance       = 1e-6;
    options.ConstraintTolerance     = 1e-6;
    options.MaxGenerations          = 30;
    options.PopulationSize          = 40;
    options.EliteCount              = ceil(0.05*options.PopulationSize);
    options.FitnessLimit            = 1e-5;
    options.InitialPopulationMatrix = []; % [0.0016, -0.0163, 0.0148]; %[0.0176, -0.0090, 0.0210];%[0.0271, 0.0048, 0.0189]; %[0.0272, 0.0048, 0.0189]; %[-0.0258, 0.0189, 0.0198];
    options.CreationFcn             = [];
    options.MaxStallGenerations     = 15;
    % options.MaxFunctionEvaluations  = 1000;
    % options.StepTolerance           = 1e-20;
    % options.FiniteDifferenceStepSize = 0.001;

    % Create problem structure to back to optimiser
    problem.fitnessfcn   =   fObj;
    problem.nvars       =   3;
    problem.options     =   options;
    problem.x0          =   p_sim_0;
    problem.Aineq       =   [];
    problem.bineq       =   [];
    problem.Aeq         =   [];
    problem.beq         =   [];
    problem.lb          =   lb;
    problem.ub          =   ub;
    problem.nonlcon     =   [];
    problem.solver      =   'ga';

    % Run ga
    [p_sim_opt, J_opt, ~, ~] = ga(problem);

    % Optimised point for DELT2
    opt_via(2,1:3) = p_sim_opt;
    clear('p_sim_opt')

    plotViaOptResults('DELT2',...
        p_sim_0,...
        opt_via(2,1:3),...
        delt2_via_downCast,...
        delt2_GP,...
        shoulder_elv,...
        osim_model,...
        data_RTSA,...
        model_MA_init.DELT2,...
        J_opt,...
        radius);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELT3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif flag_DELT3 == true
    % Search radius around init location
    radius = 0.025;
    p_sim_0 = delt3_via_loc;

    ub = p_sim_0 + radius;% delt1_via_loc + radius;%[0.05, 0.05, 0.05];
    lb = p_sim_0 - radius; %delt1_via_loc - radius;%[-0.05, -0.05, -0.05];
    lb(3) = p_sim_0(3);

    figure(101);
    hold on
    title('via point location in x-y-z space (m)')
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    scatter3(p_sim_0(1), p_sim_0(2), p_sim_0(3), 'o', 'filled','magenta')


    % Inequality constraint to keep new via point location with radius of x of
    % original via point
    %%% fCon = @(p_sim)sphere_func_con(p_sim, p_sim_0, radius);
    % Cost function to minimise moment arm differances between simulated and
    % calculated conditions
    fObj = @(p_sim)J_momentArmDist(p_sim, data_RTSA, osim_model, 'DELT3', delt3_via_downCast);

    % Set-up options
    options = optimoptions('ga', 'Display', 'iter', 'PlotFcn',{@gaplotbestf, @gaplotmaxconstr});
    options.FunctionTolerance       = 1e-6;
    options.ConstraintTolerance     = 1e-6;
    options.MaxGenerations          = 30;
    options.PopulationSize          = 40;
    options.EliteCount              = ceil(0.05*options.PopulationSize);
    options.FitnessLimit            = 1e-5;
    options.InitialPopulationMatrix = []; % [0.0016, -0.0163, 0.0148]; %[0.0176, -0.0090, 0.0210];%[0.0271, 0.0048, 0.0189]; %[0.0272, 0.0048, 0.0189]; %[-0.0258, 0.0189, 0.0198];
    options.CreationFcn             = [];
    options.MaxStallGenerations     = 15;
    % options.MaxFunctionEvaluations  = 1000;
    % options.StepTolerance           = 1e-20;
    % options.FiniteDifferenceStepSize = 0.001;

    % Create problem structure to back to optimiser
    problem.fitnessfcn   =   fObj;
    problem.nvars       =   3;
    problem.options     =   options;
    problem.x0          =   p_sim_0;
    problem.Aineq       =   [];
    problem.bineq       =   [];
    problem.Aeq         =   [];
    problem.beq         =   [];
    problem.lb          =   lb;
    problem.ub          =   ub;
    problem.nonlcon     =   [];
    problem.solver      =   'ga';

    % Run ga
    [p_sim_opt, J_opt, ~, ~] = ga(problem);

    % Optimised point for DELT1
    opt_via(3,1:3) = p_sim_opt;
    clear('p_sim_opt')

    plotViaOptResults('DELT3',...
        p_sim_0,...
        opt_via(3,1:3),...
        delt3_via_downCast,...
        delt3_GP,...
        shoulder_elv,...
        osim_model,...
        data_RTSA,...
        model_MA_init.DELT3,...
        J_opt,...
        radius);
end
%% Print out model

osim_model.print([model_file(1:end-5) '_wtf.osim'])

disp(opt_via);
end