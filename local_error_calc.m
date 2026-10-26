function [error1, error2, abs_diff, diff] = local_error(rate_func_in, sol_func, XA, t, h, BT_struct, orbit_params)
    [XB1,XB2,num_evals]=RK_step_embedded(rate_func_in,t,XA,h,BT_struct, orbit_params);
    error1 = norm(XB1 - sol_func(t + h,XA, orbit_params));
    error2 = norm(XB2 - sol_func(t + h,XA, orbit_params));
    abs_diff = norm(XB1-XB2);
    diff = norm(sol_func(t+h,XA, orbit_params) - sol_func(t,XA, orbit_params));
end