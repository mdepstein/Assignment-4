% Run your adaptive step size integrator multiple times to simulate an
% orbiting planet, across a range of desired error values. Make sure to
% record the global truncation error, the average step size, the number
% of function evaluations, and the step failure rate.
function run_explicit_variable_RK_step
    x0 = 8; y0 = 0;
    dxdt0 = 0; dydt0 = 1.5;
    V0 = [x0;y0;dxdt0;dydt0];

    h_ref = 0.05;
    t_span = [0 50];
    p = 5;
    growth_limit = 2;
    orbit_params = define_orbit_params();
    t_range = linspace(t_span(1),t_span(2),100);
    V_list = compute_planetary_motion(t_range,V0,orbit_params);

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

    h_avg_variable = zeros(size(error_desired_list));
    h_avg_fixed = zeros(size(error_desired_list));
    h_list = logspace(-2, 0, length(error_desired_list));
    num_evals_adaptive = zeros(size(error_desired_list));
    step_fail_rate_adaptive = zeros(size(error_desired_list));
    global_truncation_adaptive = zeros(size(error_desired_list));
    num_evals_fixed = zeros(size(error_desired_list));
    global_truncation_fixed = zeros(size(error_desired_list));



    % iterate through desired error for adaptive step size integrator
    for i = 1:length(error_desired_list)

        [~,X_list_adaptive,h_avg_variable(i), num_evals_adaptive(i), step_fail_rate_adaptive(i)] = ...
            explicit_RK_variable_step_integration(rate_func,t_span,V0, ...
            h_ref,DormandPrince,p,error_desired_list(i),growth_limit);
        V_actual = sol_func(t_span(end));
        V_actual = V_actual(:);
        global_truncation_adaptive(i) =  norm(X_list_adaptive(end,:)' - V_actual);
    end
    
    % iterate through step sizes for fixed step size integrator
    for j = 1:length(h_list)
        [~,X_list_fixed,h_avg_fixed(j), num_evals_fixed(j)] = explicit_RK_fixed_step_integration ...
            (rate_func,t_span,V0,h_list(j),DormandPrince, @explicit_RK_step);
        V_actual = sol_func(t_span(end));
        V_actual = V_actual(:);
        global_truncation_fixed(j) =  norm(X_list_fixed(end,:)' - V_actual);
    end

    % On the same axes, plot the global truncation error as a function of the
    % average step size for both the fixed and % adaptive step size
    % integrators. This should be a loglog plot.
    
    figure();
    % adaptive step size
    loglog(h_avg_variable, global_truncation_adaptive, 'r--');
    hold on;
    
    % fixed step size
    loglog(h_avg_fixed, global_truncation_fixed, 'b-');
    
    legend("Adaptive Step Size", "Fixed Step Size", ...
        'FontSize', 15, 'Interpreter', 'Latex');
    set(gca,'TickLabelInterpreter','latex');
    title(['Average Step Size vs. Global Truncation Error ' ...
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
        'FontSize', 15, 'Interpreter', 'Latex');
    set(gca,'TickLabelInterpreter','latex');
    title(['Number of Evaluations vs. Global Truncation Error Scaling ' ...
        '(Adapative Step Sizes)'], 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('Global Truncation Error (-)', 'Interpreter', 'latex','FontSize',15);
    xlabel('Number of Function Evaluations (-)', 'Interpreter', 'latex', 'FontSize',15);

    % Plot the failure rate as a function of the average step size, using
    % a log scale for the horizontal axis.
    figure();
    semilogx(h_avg_variable, step_fail_rate_adaptive, 'ro-');

    set(gca,'TickLabelInterpreter','latex');
    title(['Step Failure Rate vs. Average Step Size ' ...
        '(Adapative Step Sizes)'], 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('Step Failure Rate (-)', 'Interpreter', 'latex','FontSize',15);
    xlabel('Average Step Size $h$ (-)', 'Interpreter', 'latex', 'FontSize',15);

    
    % error_temp = 1e-12;
    % [t_list,X_list_adaptive,~, ~, ~] = ...
    %     explicit_RK_variable_step_integration(rate_func,t_span,V0, ...
    %     h_ref,DormandPrince,p,error_temp,growth_limit);
    % 
    % t_list2 = linspace(t_span(1),t_span(end),1001);
    % 
    % V_actual = sol_func(t_list2);
    % 
    % figure();
    % for m = 1:4
    %     subplot(4,1,m);
    %     hold on
    %     plot(t_list2,V_actual(:,m),'k','linewidth',2);
    %     plot(t_list,X_list_adaptive(:,m),'ro','markerfacecolor','r','markersize',3);
    % 
    % end


    % For a single desired error threshold, plot the position and velocity 
    % as a function of time (position vs. time and velocity vs. time 
    % should be on separate axes). Do this with both lines and dots:
    
    error_desired = error_desired_list(15);
    error_desired_str = sprintf('%.1e', error_desired);

    [t_list3,X_list_single_error_thres,~, ~, ~] = explicit_RK_variable_step_integration ...
        (rate_func,t_span,V0,h_ref,DormandPrince,p,error_desired,growth_limit);

    figure();
    subplot(2,1,1);
    plot(t_list3, X_list_single_error_thres(:,1), 'ro-');
    hold on;
    plot(t_list3, X_list_single_error_thres(:,2), 'bo-');

    set(gca,'TickLabelInterpreter','latex');
    title(['Position vs. Time, error threshold = ' error_desired_str], ...
        'Interpreter', 'latex', 'FontSize', 15);
    ylabel('Position (-)', 'Interpreter', 'latex','FontSize',15);
    xlabel('Time (-)', 'Interpreter', 'latex', 'FontSize',15);
    legend('$x$-position', '$y$-position', ...
        'Interpreter', 'latex', 'Location', 'best');
    hold off;

    subplot(2,1,2);
    plot(t_list3, X_list_single_error_thres(:,3), 'ro-');
    hold on;
    plot(t_list3, X_list_single_error_thres(:,4), 'bo-');
    
    set(gca,'TickLabelInterpreter','latex');
    title(['Velocity as a Function of Time, error threshold = ' error_desired_str], ...
        'Interpreter', 'latex', 'FontSize', 15);
    ylabel('Velocity (-)', 'Interpreter', 'latex','FontSize',15);
    xlabel('Time (-)', 'Interpreter', 'latex', 'FontSize',15);
    legend('$x$-velocity', '$y$-velocity', ...
        'Interpreter', 'latex', 'Location', 'best');
    hold off;


    % For a single desired error threshold, generate a scatter plot of the 
    % step size as a function of the distance between the planet and the
    % sun, r = x2 +y2. Note: you can extract the list of step sizes by just 
    % taking the difference between adjacent times
    
    sun_dist = zeros(length(X_list_single_error_thres)-1, 1);

    for k = 1:length(X_list_single_error_thres)-1
        x = X_list_single_error_thres(k,1);
        y = X_list_single_error_thres(k,2);
        sun_dist(k) = sqrt(x^2+y^2);
    end
    
    step_sizes = diff(t_list3);
    
    figure();
    scatter(sun_dist, step_sizes');

    set(gca,'TickLabelInterpreter','latex');
    title(['Adaptive Step Size vs. Distance from the Sun, ' ...
        'error threshold = ' error_desired_str], ...
        'Interpreter', 'latex', 'FontSize', 15);
    ylabel('Step Size $h$ (-)', 'Interpreter', 'latex','FontSize',15);
    xlabel('Distance from the Sun (-)', 'Interpreter', 'latex', 'FontSize',15);

end