function run_method()
figure()
orbit_params = define_orbit_params();
BT_struct = butcher_tableau_struct();
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
h_ref = .3;
X0 = V_list(1,1:4)
color = ["blue","red","green","magenta","orange"];
BT_struct_names = ["DormandPrince", "Fehlberg", "FehlbergRK1", "HeunEuler", "Bogacki"];
for i = 1:3
[t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(@gravity_rate_func,t_span,X0,h_ref,BT_struct, i)
plot(t_list, X_listlist, 'Marker', 'o','Color',color(i),'MarkerSize',2,'MarkerFaceColor',color(i));
end
l = legend("Analytical Solution",BT_struct_names(1),BT_struct_names(2),BT_struct_names(3),'FontSize',10 );



end