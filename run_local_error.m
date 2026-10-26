function run_local_error(index)
    %close all;
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    orbit_params = define_orbit_params();

    rate_func = @(tin,Vin) gravity_rate_func(tin, Vin, orbit_params);
    

    x0 = 8;
    y0 = .2;
    dxdt0 = -.3;
    dydt0 = 1.5;
    V0 = [x0;y0;dxdt0;dydt0];
    
    BT_struct_var = butcher_tableau_struct_var_step();
    BT_struct_var_names = ["DormandPrince", "Fehlberg", "FehlbergRK1", "HeunEuler", "Bogacki"];
    
    BT_struct = BT_struct_var(index);
    t_ref = 7;
    sol_func = @(tin) compute_planetary_motion(tin, V0, orbit_params);
    VA = sol_func(t_ref);
    
    %XA = solution01(t,0,orbit_params); %0 is just an arbituary value not needed for calculating solution
    hs = logspace(-1,0,100);
    local_errors_1 = zeros(size(hs));
    local_errors_2 = zeros(size(hs));
    approx_sol_diff = zeros(size(hs));
    solution_diff = zeros(size(hs));
    
    for n = 1:length(hs)
        h = hs(n);
        [error1, error2, abs_diff, diff] = local_error_calc(rate_func, sol_func, VA, t_ref, h, BT_struct);
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

    loglog(hs,exp(log(hs)*p_XB1(1)+p_XB1(2)),'r','linewidth',2);
    loglog(hs,exp(log(hs)*p_XB2(1)+p_XB2(2)),'r','linewidth',2);
    % loglog(hs,exp(log(hs)*p_XB1(1)+p_XB1(2)),'r','linewidth',2);


    myfig = gcf();
    exportgraphics(myfig , 'local_err_explicit_01.png', 'Resolution', 600);
    




end
