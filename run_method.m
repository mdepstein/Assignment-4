function run_method()
figure()
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
axis equal; axis square;
axis([-20,20,-20,20])
hold on
plot(0,0,'ro','markerfacecolor','r','markersize',5);
plot(V_list(:,1),V_list(:,2),'k');
h_ref = .1;
X0 = V_list(1,1:4)
color = ["blue","red","green","magenta","orange"];
BT_struct_var_names = ["DormandPrince", "Fehlberg", "FehlbergRK1", "HeunEuler", "Bogacki"];
BT_struct_explicit_names = ["Forward", "Midpoint", "Kutta"];

for i = 1:3
[t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(@gravity_rate_func,t_span,X0,h_ref,BT_struct_explicit(i));
size(X_list)
plot(X_list(1,:), X_list(2,:), 'Marker', 'o','Color',color(i),'MarkerSize',2,'MarkerFaceColor',color(i));
end
l = legend("Sun", "Analytical Solution",BT_struct_explicit_names(1),BT_struct_explicit_names(2),BT_struct_explicit_names(3),'FontSize',10 );



end