function run_method()
figure()
 set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
orbit_params = define_orbit_params();
BT_struct_var = butcher_tableau_struct_var_step();
BT_struct_explicit = butcher_tableau_struct_explicit();
x0 = 8;
y0 = 0;
dxdt0 = 0;
dydt0 = 1.5;
V0 = [x0;y0;dxdt0;dydt0];
t_span = [0,50];
t_range = linspace(t_span(1),t_span(2),100);
V_list = compute_planetary_motion(t_range,V0,orbit_params);
h_ref = .05;

%plot(0,0,'ro','markerfacecolor','r','markersize',5);

X0 = V_list(1,1:4)
color = ["blue","red","green","magenta","orange"];
BT_struct_var_names = ["DormandPrince", "Fehlberg", "FehlbergRK1", "HeunEuler", "Bogacki"];
BT_struct_explicit_names = ["Ralston's", "SSPRK3", "RK4"];
t = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(1)
plot(t_range, V_list(:,1),'k');
hold on
nexttile(2)
plot(t_range,V_list(:,2),'k');
hold on

wrap_gravity_rate_func = @(t,XA) gravity_rate_func(t,XA,orbit_params);
for i = 1:3
    nexttile(1)
    [t_list,X_list,h_avg, num_evals] = ...
        explicit_RK_fixed_step_integration(wrap_gravity_rate_func,t_span,X0,h_ref, BT_struct_explicit(i), @explicit_RK_step);
    plot(t_list, X_list(1,:), 'Marker', 'o','Color',color(i),'MarkerSize',2,'MarkerFaceColor',color(i));
    xlabel("Time (-)",'FontSize',15)
        ylabel("$x(t) (-)$",'FontSize',15)
        axis([t_span(1) t_span(2) -5 9])
        ax = gca; 
        ax.FontSize = 14; 

    hold on
    nexttile(2)
    plot(t_list, X_list(2,:),'Marker', 'o','Color',color(i), 'MarkerSize',2,'MarkerFaceColor',color(i));

    xlabel("Time (-)",'FontSize',15)
        ylabel("$y(t) (-)$",'FontSize',15)
        axis([t_span(1),t_span(2) -5 9])
        ax = gca; 
        ax.FontSize = 14; 
    hold on
end


title(t, 'Comparison of Runge Kutta Explicit Methods Integration Plot', 'FontSize', 15, 'Interpreter', 'latex')    

l = legend("Analytical Solution", ...
    BT_struct_explicit_names(1), ...
    BT_struct_explicit_names(2), ...
    BT_struct_explicit_names(3), ...
    'FontSize', 10, 'Interpreter', 'Latex');

set(gca,'TickLabelInterpreter','latex');
l.Layout.Tile = 'north'; 
    myfig = gcf();
    exportgraphics(myfig , 'explicit_integration.png', 'Resolution', 600);


end