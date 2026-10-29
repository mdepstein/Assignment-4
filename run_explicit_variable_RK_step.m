% Run your adaptive step size integrator multiple times to simulate an 
% orbiting planet, across a range of desired error values. Make sure to
% record the global truncation error, the average step size, the number 
% of function evaluations, and the step failure rate.
function run_explicit_variable_RK_step
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    V0 = [x0;y0;dxdt0;dydt0];
    h_ref = 0.3;
    tspan = [0 20];
    p = 5;
    growth_limit = 2;
    orbit_params = define_orbit_params();
    rate_func = @(tin,Vin) gravity_rate_func(tin, Vin, orbit_params);
    sol_func = @(tin) compute_planetary_motion(tin, V0, orbit_params);

    
    error_desired_list = logspace(-5, -1, 50);
    DormandPrince = struct();
    DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    DormandPrince.B = [35/384, 0, 500/1113, 125/192,-2187/6784, 11/84, 0;...
        5179/57600, 0, 7571/16695, 393/640,-92097/339200, 187/2100, 1/40];
    DormandPrince.A = [0,0,0,0,0,0,0;
        1/5, 0, 0, 0,0,0,0;...
        3/40, 9/40, 0, 0, 0, 0,0;...
        44/45,-56/15, 32/9, 0, 0, 0,0;...
        19372/6561,-25360/2187, 64448/6561,-212/729, 0, 0,0;...
        9017/3168,-355/33, 46732/5247, 49/176,-5103/18656, 0,0;...
        35/384, 0, 500/1113, 125/192,-2187/6784, 11/84,0];

    h_avg = zeros(size(error_desired_list));
    h_list = logspace(-2, 0, length(error_desired_list));
    num_evals_adaptive = zeros(size(error_desired_list));
    step_fail_rate_adaptive = zeros(size(error_desired_list));
    global_truncation_adaptive = zeros(size(error_desired_list));
    num_evals_fixed = zeros(size(error_desired_list));
    global_truncation_fixed = zeros(size(error_desired_list));

    % iterate through desired error for adaptive step size integrator
    for i = 1:length(error_desired_list)

        [~,X_list_adaptive,h_avg(i), num_evals_adaptive(i), step_fail_rate_adaptive(i)] = ...
            explicit_RK_variable_step_integration(rate_func,tspan,V0, ...
            h_ref,DormandPrince,p,error_desired_list(i),growth_limit);
        V_actual = sol_func(tspan(end));
        global_truncation_adaptive(i) =  norm(X_list_adaptive(end,:) - V_actual);
    end
    
    % iterate through step sizes for fixed step size integrator
    for j = 1:length(h_list)
        [~,X_list_fixed,~, num_evals_fixed(j)] = explicit_RK_fixed_step_integration ...
            (@gravity_rate_func,tspan,V0,h_list(j),DormandPrince, @explicit_RK_step, orbit_params);
        V_actual = sol_func(tspan(end));
        global_truncation_fixed(j) =  norm(X_list_fixed(end,:) - V_actual);
    end

    % On the same axes, plot the global truncation error as a function of the 
    % average step size for both the fixed and % adaptive step size 
    % integrators. This should be a loglog plot.
    
    figure();
    % adaptive step size
    loglog(h_avg, global_truncation_adaptive, 'r--');
    hold on;
    
    % fixed step size
    loglog(h_list, global_truncation_fixed, 'b-');

    legend("Adaptive Step Size", "Fixed Step Size", ...
        'FontSize', 10, 'Interpreter', 'Latex');
    set(gca,'TickLabelInterpreter','latex');
    title(['Average Step Size vs. Global Truncation Error' ...
        '(Adapative Step Sizes)'], 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('Global Truncation Error (-)', 'Interpreter', 'latex','FontSize',15);
    xlabel('Average Step Size $h$ (-)', 'Interpreter', 'latex', 'FontSize',15);


    %  On the same axes, plot the global truncation error as a function of 
    % the number of function evaluations for both the fixed and adaptive 
    % step size integrators.
    
    figure();
    % adaptive step size
    loglog(num_evals_adaptive, global_truncation_adaptive, 'r--');
    hold on;

    % fixed step size
    loglog(num_evals_fixed, global_truncation_fixed, 'b-');

    legend("Adaptive Step Size", "Fixed Step Size", ...
        'FontSize', 10, 'Interpreter', 'Latex');
    set(gca,'TickLabelInterpreter','latex');
    title(['Average Step Size vs. Global Truncation Error Scaling' ...
        '(Adapative Step Sizes)'], 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('Global Truncation Error (-)', 'Interpreter', 'latex','FontSize',15);
    xlabel('Number of Function Evaluations (-)', 'Interpreter', 'latex', 'FontSize',15);

    % Plot the failure rate as a function of the average step size, using 
    % a log scale for the horizontal axis.
    figure();
    loglog(h_avg, step_fail_rate_adaptive, 'r--');

    set(gca,'TickLabelInterpreter','latex');
    title(['Average Step Size vs. Failure Rate' ...
        '(Adapative Step Sizes)'], 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('Failure Rate (-)', 'Interpreter', 'latex','FontSize',15);
    xlabel('Average Step Size $h$ (-)', 'Interpreter', 'latex', 'FontSize',15);

end