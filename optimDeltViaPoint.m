function optimDeltViaPoint(model_file, flag_DELT1, flag_DELT2, flag_DELT3)
%optimDeltViaPoint Defines DELT1, DELT2 and DELT3 via poit locations
%   Detailed explanation goes here

%% Set-up
% GA setup
function_tolerance       = 1e-6;
constraint_colerance     = 1e-6;
max_generations          = 75;
population_size          = 120;
elite_count              = ceil(0.05*population_size);
fitness_limit            = 1e-5;
initial_populationMatrix = []; % DELT3 insertion [0.0016, -0.0163, 0.0148]; %[0.0176, -0.0090, 0.0210];%[0.0271, 0.0048, 0.0189]; %[0.0272, 0.0048, 0.0189]; %[-0.0258, 0.0189, 0.0198];
creation_fcn             = [];
max_stallGenerations     = 20;

% Flags
flag_delt1ViaDelete = false;
flag_delt2ViaDelete = false;
flag_delt3ViaDelete = false;

% Get model
import org.opensim.modeling.*

osim_model = Model(model_file);
init_state = osim_model.initSystem();

% Ackland et al (2010) RTSA MA data

data_adb_RTSA.angles        = [2.5, 30, 60, 90, 120];
data_adb_RTSA.DELT1_mean    = [15.6, 25.2,32.5, 35.8, 33.3]*0.001;
data_adb_RTSA.DELT1_sd      = [2.3, 2.5, 1.8, 3.4, 2.7]*0.001;
data_adb_RTSA.DELT2_mean    = [30.2, 33.9, 42.2, 46.2, 39.8]*0.001;
data_adb_RTSA.DELT2_sd      = [6.6, 5.5, 5.5, 4.2, 5.8]*0.001;
data_adb_RTSA.DELT3_mean    = [1.3, 3.5, 7.3, 11.4, 14.1]*0.001;
data_adb_RTSA.DELT3_sd      = [1.4, 1.0, 1.7, 2.5, 3.6]*0.001;

data_flx_RTSA.angles        = [2.5, 30, 60, 90, 120];
data_flx_RTSA.DELT1_mean    = [25.9, 34.2, 35.8, 35.7, 32.7]*0.001;
data_flx_RTSA.DELT1_sd      = [5.1, 2.9, 4.4, 5.3, 5.6]*0.001;
data_flx_RTSA.DELT2_mean    = [14.2, 18.6, 20.1, 22.8, 27.0]*0.001;
data_flx_RTSA.DELT2_sd      = [6.4, 7.1, 4.9, 6.6, 6.1]*0.001;
data_flx_RTSA.DELT3_mean    = [-14.6, -17.6, -16.3, -13.8, -13.1]*0.001;
data_flx_RTSA.DELT3_sd      = [5.5, 7.0, 5.1, 2.8, 2.4]*0.001;

%% Handle model

% Get muscle handles
delt1 = osim_model.getMuscles.get('DELT1');
delt2 = osim_model.getMuscles.get('DELT2');
delt3 = osim_model.getMuscles.get('DELT3');

% Get coordinate handle
shoulder_elv = osim_model.getCoordinateSet().get('shoulder_elv');
elv_angle    = osim_model.getCoordinateSet().get('elv_angle');

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
for i_angle = 1:length(data_adb_RTSA.angles)

    osim_model.updCoordinateSet().get('elv_angle').setValue(init_state, 0, true);
    osim_model.updCoordinateSet().get('shoulder_elv').setValue(init_state, deg2rad(data_adb_RTSA.angles(i_angle)));
    osim_model.realizePosition(init_state);

    model_abd_MA_init.DELT1(i_angle) = delt1_GP.computeMomentArm(init_state,shoulder_elv);

    osim_model.updCoordinateSet().get('elv_angle').setValue(init_state, deg2rad(70), true);
    osim_model.updCoordinateSet().get('shoulder_elv').setValue(init_state, deg2rad(data_flx_RTSA.angles(i_angle)));
    osim_model.realizePosition(init_state);

    model_flx_MA_init.DELT1(i_angle) = delt1_GP.computeMomentArm(init_state,elv_angle);

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
        
        point_count_3 = i_point;

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
for i_angle = 1:length(data_adb_RTSA.angles)

    osim_model.updCoordinateSet().get('elv_angle').setValue(init_state, 0, true);
    osim_model.updCoordinateSet().get('shoulder_elv').setValue(init_state, deg2rad(data_adb_RTSA.angles(i_angle)));
    osim_model.realizePosition(init_state);

    model_abd_MA_init.DELT2(i_angle) = delt2_GP.computeMomentArm(init_state,shoulder_elv);

    osim_model.updCoordinateSet().get('elv_angle').setValue(init_state, deg2rad(70), true);
    osim_model.updCoordinateSet().get('shoulder_elv').setValue(init_state, deg2rad(data_flx_RTSA.angles(i_angle)));
    osim_model.realizePosition(init_state);

    model_flx_MA_init.DELT2(i_angle) = delt2_GP.computeMomentArm(init_state,elv_angle);

end

% Delete via point
if flag_delt2ViaDelete == true
   delt2_GP.deletePathPoint(init_state, point_count_3);
   osim_model.initSystem();
   osim_model.finalizeConnections();
   osim_model.initSystem();
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

        point_count_3 = i_point;

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
for i_angle = 1:length(data_adb_RTSA.angles)

    osim_model.updCoordinateSet().get('elv_angle').setValue(init_state, 0, true);
    osim_model.updCoordinateSet().get('shoulder_elv').setValue(init_state, deg2rad(data_adb_RTSA.angles(i_angle)));
    osim_model.realizePosition(init_state);

    model_abd_MA_init.DELT3(i_angle) = delt3_GP.computeMomentArm(init_state,shoulder_elv);

    osim_model.updCoordinateSet().get('elv_angle').setValue(init_state, deg2rad(70), true);
    osim_model.updCoordinateSet().get('shoulder_elv').setValue(init_state, deg2rad(data_flx_RTSA.angles(i_angle)));
    osim_model.realizePosition(init_state);

    model_flx_MA_init.DELT3(i_angle) = delt3_GP.computeMomentArm(init_state,elv_angle);

end

% Delete via point
if flag_delt3ViaDelete == true
   delt3_GP.deletePathPoint(init_state, point_count_3);
   osim_model.initSystem();
   osim_model.finalizeConnections();
   osim_model.initSystem();
end
%% Optimisation - MA calculation after via point and joint modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELT1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_DELT1 == true
    % Search radius around init location
    radius = 0.035; %0.025
    p_sim_0 = delt1_via_loc;

    ub = p_sim_0 + radius;% delt1_via_loc + radius;%[0.05, 0.05, 0.05];
    lb = p_sim_0 - radius; %delt1_via_loc - radius;%[-0.05, -0.05, -0.05];
%     lb(3) = p_sim_0(3);

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
    fObj = @(p_sim)J_momentArmDist(p_sim, data_adb_RTSA, data_flx_RTSA, osim_model, 'DELT1', delt1_via_downCast, flag_delt1ViaDelete);

    % Set-up options
    options = optimoptions('ga', 'Display', 'iter', 'PlotFcn',{@gaplotbestf, @gaplotmaxconstr});
    options.FunctionTolerance       = function_tolerance;
    options.ConstraintTolerance     = constraint_colerance;
    options.MaxGenerations          = max_generations;
    options.PopulationSize          = population_size;
    options.EliteCount              = elite_count;
    options.FitnessLimit            = fitness_limit;
    options.InitialPopulationMatrix = initial_populationMatrix; 
    options.CreationFcn             = creation_fcn;
    options.MaxStallGenerations     = max_stallGenerations;
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
        data_adb_RTSA,...
        model_abd_MA_init.DELT1,...
        J_opt,...
        radius,...
        flag_delt1ViaDelete);

    plotViaOptResults('DELT1',...
        p_sim_0,...
        opt_via(1,1:3),...
        delt1_via_downCast,...
        delt1_GP,...
        elv_angle,...
        osim_model,...
        data_flx_RTSA,...
        model_flx_MA_init.DELT1,...
        J_opt,...
        radius,...
        flag_delt1ViaDelete);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELT2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_DELT2 == true
    % Search radius around init location
    radius = 0.035;
    p_sim_0 = delt2_via_loc;

    ub = p_sim_0 + radius;% delt1_via_loc + radius;%[0.05, 0.05, 0.05];
    lb = p_sim_0 - radius; %delt1_via_loc - radius;%[-0.05, -0.05, -0.05];
%     ub(1) = p_sim_0(1);
%     ub(2) = p_sim_0(2);
%     lb(3) = p_sim_0(3);

    figure(101);
    hold on
    title('via point location in x-y-z space (m)')
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    scatter3(p_sim_0(1), p_sim_0(2), p_sim_0(3), 'o', 'filled','magenta')



    if flag_delt2ViaDelete == true
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Optimised point for DELT2 - THIS IS DELT2 ORIGIN (needs deletion
        % but best workaround for best [physiological] result)
        opt_via(2,1:3) = [0, 0, 0];

        error_MA_DELT2 = J_momentArmDist(opt_via(2,1:3), data_adb_RTSA, data_flx_RTSA, osim_model, 'DELT2', delt2_via_downCast, flag_delt2ViaDelete);

        plotViaOptResults('DELT2',...
            p_sim_0,...
            opt_via(2,1:3),...
            delt2_via_downCast,...
            delt2_GP,...
            shoulder_elv,...
            osim_model,...
            data_adb_RTSA,...
            model_abd_MA_init.DELT2,...
            error_MA_DELT2,...
            radius,...
            flag_delt2ViaDelete);

        plotViaOptResults('DELT2',...
            p_sim_0,...
            opt_via(2,1:3),...
            delt2_via_downCast,...
            delt2_GP,...
            elv_angle,...
            osim_model,...
            data_flx_RTSA,...
            model_flx_MA_init.DELT2,...
            error_MA_DELT2,...
            radius,...
            flag_delt2ViaDelete);
    else
        % Inequality constraint to keep new via point location with radius of x of
        % original via point
        %%% fCon = @(p_sim)sphere_func_con(p_sim, p_sim_0, radius);
        % Cost function to minimise moment arm differances between simulated and
        % calculated conditions
        fObj = @(p_sim)J_momentArmDist(p_sim, data_adb_RTSA, data_flx_RTSA, osim_model, 'DELT2', delt2_via_downCast, flag_delt2ViaDelete);

        % Set-up options
        options = optimoptions('ga', 'Display', 'iter', 'PlotFcn',{@gaplotbestf, @gaplotmaxconstr});
        options.FunctionTolerance       = function_tolerance;
        options.ConstraintTolerance     = constraint_colerance;
        options.MaxGenerations          = max_generations;
        options.PopulationSize          = population_size;
        options.EliteCount              = elite_count;
        options.FitnessLimit            = fitness_limit;
        options.InitialPopulationMatrix = initial_populationMatrix;
        options.CreationFcn             = creation_fcn;
        options.MaxStallGenerations     = max_stallGenerations;
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
        opt_via(2,1:3) = p_sim_opt;
        clear('p_sim_opt')

        plotViaOptResults('DELT2',...
            p_sim_0,...
            opt_via(2,1:3),...
            delt2_via_downCast,...
            delt2_GP,...
            shoulder_elv,...
            osim_model,...
            data_adb_RTSA,...
            model_abd_MA_init.DELT2,...
            J_opt,...
            radius,...
            flag_delt2ViaDelete);

        plotViaOptResults('DELT2',...
            p_sim_0,...
            opt_via(2,1:3),...
            delt2_via_downCast,...
            delt2_GP,...
            elv_angle,...
            osim_model,...
            data_flx_RTSA,...
            model_flx_MA_init.DELT2,...
            J_opt,...
            radius,...
            flag_delt2ViaDelete);
    end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELT3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_DELT3 == true
    % Search radius around init location
    radius = 0.035;
    p_sim_0 = delt3_via_loc;

    ub = p_sim_0 + radius;% delt1_via_loc + radius;%[0.05, 0.05, 0.05];
    lb = p_sim_0 - radius; %delt1_via_loc - radius;%[-0.05, -0.05, -0.05];
%     ub(1) = p_sim_0(1);
    %     lb(3) = p_sim_0(3);

    figure(101);
    hold on
    title('via point location in x-y-z space (m)')
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    scatter3(p_sim_0(1), p_sim_0(2), p_sim_0(3), 'o', 'filled','magenta')





    if flag_delt3ViaDelete == true
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Optimised point for DELT3 - THIS IS DELT3 ORIGIN (needs deletion
        % but best workaround for best [physiological] result)
        opt_via(3,1:3) = [0, 0, 0];

        error_MA_DELT3 = J_momentArmDist(opt_via(3,1:3), data_adb_RTSA, data_flx_RTSA, osim_model, 'DELT3', delt3_via_downCast, flag_delt3ViaDelete);

        plotViaOptResults('DELT3',...
            p_sim_0,...
            opt_via(3,1:3),...
            delt3_via_downCast,...
            delt3_GP,...
            shoulder_elv,...
            osim_model,...
            data_adb_RTSA,...
            model_abd_MA_init.DELT3,...
            error_MA_DELT3,...
            radius,...
            flag_delt3ViaDelete);

        plotViaOptResults('DELT3',...
            p_sim_0,...
            opt_via(3,1:3),...
            delt3_via_downCast,...
            delt3_GP,...
            elv_angle,...
            osim_model,...
            data_flx_RTSA,...
            model_flx_MA_init.DELT3,...
            error_MA_DELT3,...
            radius,...
            flag_delt3ViaDelete);
    else

        % Inequality constraint to keep new via point location with radius of x of
        % original via point
        %%% fCon = @(p_sim)sphere_func_con(p_sim, p_sim_0, radius);
        % Cost function to minimise moment arm differances between simulated and
        % calculated conditions
        fObj = @(p_sim)J_momentArmDist(p_sim, data_adb_RTSA, data_flx_RTSA, osim_model, 'DELT3', delt3_via_downCast, flag_delt2ViaDelete);

        % Set-up options
        options = optimoptions('ga', 'Display', 'iter', 'PlotFcn',{@gaplotbestf, @gaplotmaxconstr});
        options.FunctionTolerance       = function_tolerance;
        options.ConstraintTolerance     = constraint_colerance;
        options.MaxGenerations          = max_generations;
        options.PopulationSize          = population_size;
        options.EliteCount              = elite_count;
        options.FitnessLimit            = fitness_limit;
        options.InitialPopulationMatrix = initial_populationMatrix;
        options.CreationFcn             = creation_fcn;
        options.MaxStallGenerations     = max_stallGenerations;
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
            data_adb_RTSA,...
            model_abd_MA_init.DELT3,...
            J_opt,...
            radius,...
            flag_delt2ViaDelete);

        plotViaOptResults('DELT3',...
            p_sim_0,...
            opt_via(3,1:3),...
            delt3_via_downCast,...
            delt3_GP,...
            elv_angle,...
            osim_model,...
            data_flx_RTSA,...
            model_flx_MA_init.DELT3,...
            J_opt,...
            radius,...
            flag_delt3ViaDelete);
    end

end
%% Print out model and save figures
osim_model.print(model_file);
% osim_model.print([model_file(1:end-5) '_wtf.osim']);

saveas(figure(2),['C:\Users\lab\Pavlos\Research\UVic\Project_RTSA_config_eval\OpenSim\Out\Verification_Validation\Moment_arms\DELT1_shoulder_elv_10-10_' model_file(end-15:end-5) '.tif'],'tif')
saveas(figure(3),['C:\Users\lab\Pavlos\Research\UVic\Project_RTSA_config_eval\OpenSim\Out\Verification_Validation\Moment_arms\DELT1_elv_angle_10-10_' model_file(end-15:end-5) '.tif'],'tif')
saveas(figure(4),['C:\Users\lab\Pavlos\Research\UVic\Project_RTSA_config_eval\OpenSim\Out\Verification_Validation\Moment_arms\DELT2_shoulder_elv_10-10_' model_file(end-15:end-5) '.tif'],'tif')
saveas(figure(5),['C:\Users\lab\Pavlos\Research\UVic\Project_RTSA_config_eval\OpenSim\Out\Verification_Validation\Moment_arms\DELT2_elv_angle_10-10_' model_file(end-15:end-5) '.tif'],'tif')
saveas(figure(6),['C:\Users\lab\Pavlos\Research\UVic\Project_RTSA_config_eval\OpenSim\Out\Verification_Validation\Moment_arms\DELT3_shoulder_elv_10-10_' model_file(end-15:end-5) '.tif'],'tif')
saveas(figure(7),['C:\Users\lab\Pavlos\Research\UVic\Project_RTSA_config_eval\OpenSim\Out\Verification_Validation\Moment_arms\DELT3_elv_angle_10-10_' model_file(end-15:end-5) '.tif'],'tif')
disp(opt_via);
end