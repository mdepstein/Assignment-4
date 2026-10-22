function global_truncation
    orbit_params = define_orbit_params();
    BT_struct_var = butcher_tableau_struct_var_step();
    BT_struct_explicit = butcher_tableau_struct_explicit();
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    V0 = [x0;y0;dxdt0;dydt0];
    t_span = [0,30];
    t_range = linspace(t_span(1),t_span(2),100);
    V_list = compute_planetary_motion(t_range,V0,orbit_params);
    h_list = logspace(-3, -2, 20);
    global_truncation = zeros(3, length(h_list));

    % X0 = V_list(1,1:4);
    X0=V0
    color = ["blue","red","green"];
    BT_struct_var_names = ["DormandPrince", "Fehlberg", "FehlbergRK1", "HeunEuler", "Bogacki"];
    BT_struct_explicit_names = ["Forward Euler", "Midpoint Method", "Runga Kutta"];

    for i = 1:3
        for j = 1:length(h_list)
            h_ref = h_list(j);

            [t_list,X_list,h_avg, num_evals] = ...
                explicit_RK_fixed_step_integration_global( ...
                @gravity_rate_func,t_span,X0,h_ref, ...
                BT_struct_explicit(i));
            V_actual = compute_planetary_motion(t_list(end),V0,orbit_params);

            % size(X_list)
            % size(V_actual)
            global_truncation(i, j) = norm(X_list(:, end) - V_actual, 2);
        end
    end

    figure;
    hold on;

    for k = 1:3
        loglog(h_list, global_truncation(k, :), ...
        '.', 'MarkerSize', 15, 'Color', color(k));
        fit = polyfit(log(h_list), log(global_truncation(k, :)), 1);
        fit_error = exp(polyval(fit, log(h_list)));
        loglog(h_list, fit_error, '--', 'color', color(k));
        p(:, k) = [BT_struct_explicit_names(k), (fit(1))];
    end
    p


    l = legend("Forward Euler", "Forward Euler Fit Line", ...
    "Midpoint Method", "Midpoint Method Fit Line", "Runga Kutta", ...
    "Runga Kutta Fit Line", 'FontSize', 10, ...
    'Interpreter', 'Latex', 'Location', 'best');

    set(gca,'TickLabelInterpreter','latex');
    title('Global Truncation Error vs. Step Size', ...
        'Interpreter', 'latex', 'FontSize', 15);
    xlabel('Step Size $h_{ref}$ (-)', 'Interpreter', 'latex');
    ylabel('Global Truncation (-)', 'Interpreter', 'latex');
end