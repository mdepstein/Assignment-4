function run_local_error()
    %close all;
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    orbit_params = define_orbit_params();
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    V0 = [x0;y0;dxdt0;dydt0];
    
    BT_struct_var = butcher_tableau_struct_var_step();
    BT_struct_var_names = ["DormandPrince", "Fehlberg", "FehlbergRK1", "HeunEuler", "Bogacki"];
    index= 2;
    BT_struct = BT_struct_var(index);
    t = 7;
    VA = compute_planetary_motion(t,V0,orbit_params);
    VA =
    %XA = solution01(t,0,orbit_params); %0 is just an arbituary value not needed for calculating solution
    hs = logspace(-4.5,.5,100);
    local_errors_1 = zeros(size(hs));
    local_errors_2 = zeros(size(hs));
    approx_sol_diff = zeros(size(hs));
    solution_diff = zeros(size(hs));
    
    for n = 1:length(hs)
        h = hs(n);
        [error1, error2, abs_diff, diff] = local_error_calc(@gravity_rate_func, @compute_planetary_motion, VA, t, h, BT_struct, orbit_params);
        %[error1, error2, abs_diff, diff] = local_error_calc(@rate_func01, @solution01, XA, t, h, BT_struct, orbit_params);
        local_errors_1(n) = error1;
        local_errors_2(n) = error2;
        approx_sol_diff(n) = abs_diff;
        solution_diff(n) = diff;

        
    end

    figure();
    loglog(hs, local_errors_1, "b.", MarkerSize=13); hold on
    loglog(hs, local_errors_2, "c.", MarkerSize=10);

    loglog(hs, approx_sol_diff, "k.", MarkerSize=10); hold on
    loglog(hs, solution_diff, "m.", MarkerSize=10);

    title("Local Truncation Error for Embedded RK Method",'FontSize',13);
    subtitle(BT_struct_var_names(index));
    xlabel("$h_{avg} (-)$ ",'FontSize',13);
    ylabel("Local Truncation error $(-)$",'FontSize',13)
    legend("XB1 Local Error", "XB2 Local Error", "|XB1-XB2|", "Ground Truth", 'FontSize',11, "Location", "southeast")
    
    p_XB1 = polyfit(log(hs), log(local_errors_1), 1)
    p_XB2 = polyfit(log(hs), log(local_errors_2), 1)
    p_diff = polyfit(log(hs), log(approx_sol_diff), 1)

    myfig = gcf();
    exportgraphics(myfig , 'local_err_explicit_01.png', 'Resolution', 600);
    




end
