function run_local_error()
    close all;
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');

    t = 0.492;
    compute_planetary_motion(t_list,V0,orbit_params)    
    hs = logspace(-4.5,-1,100);
    local_errors_e = zeros(size(hs));
    local_errors_m = zeros(size(hs));
    local_errors_b = zeros(size(hs));
    local_errors_i = zeros(size(hs));

    diffs_e = zeros(size(hs));
    diffs_m = zeros(size(hs));
    
    for n = 1:length(hs)
        h = hs(n);
        [l_ee, diffe] = local_error(@rate_func01, @solution01, XA, t, h, "euler");
        [l_em, diffm] = local_error(@rate_func01, @solution01, XA, t, h, "midpoint");
        local_errors_e(n) = l_ee;
        local_errors_m(n) = l_em;
        diffs_e(n) = diffe;
        diffs_m(n) = diffm;

        [l_eb, diffb] = local_error(@rate_func01, @solution01, XA, t, h, "backward_euler");
        [l_ei, diffi] = local_error(@rate_func01, @solution01, XA, t, h, "implicit_midpoint");
        local_errors_b(n) = l_eb;
        local_errors_i(n) = l_ei;
        
    end
    
    figure();
    loglog(hs, local_errors_e, "b.", MarkerSize=10); hold on
    loglog(hs, diffs_e, "c*", MarkerSize=10);
    coeffs_e = polyfit(log(hs), log(local_errors_e), 1);
    localfit_e = exp(polyval(coeffs_e, log(hs)));
    loglog(hs, localfit_e, "r", LineWidth=1);

    loglog(hs, local_errors_m, "k.", MarkerSize=10); hold on
    loglog(hs, diffs_m, "m.", MarkerSize=10);
    coeffs_m = polyfit(log(hs), log(local_errors_m), 1);
    localfit_m = exp(polyval(coeffs_m, log(hs)));
    loglog(hs, localfit_m, "g", LineWidth=1);

    % loglog(hs, local_errors_b, "m.", MarkerSize=10); hold on
    % 
    % loglog(hs, local_errors_i, "c.", MarkerSize=10); hold on

    title("Local Truncation Error for Explicit Methods",'FontSize',13);
    xlabel("$h_{avg} (-)$ ",'FontSize',13);
    ylabel("Local Truncation error $(-)$",'FontSize',13)
    legend("Euler method error", "Euler differences", "Euler fit line", "Midpoint method error", "Midpoint differences", "Midpoint fit line", 'FontSize',11, "Location", "southeast")

    myfig = gcf();
    exportgraphics(myfig , 'local_err_explicit_01.png', 'Resolution', 600);