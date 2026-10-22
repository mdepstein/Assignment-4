function [error, diff] = local_error(rate_func_in, sol_func, XA, t, h, method)
    Xn1 = method(rate_func_in, t, XA, h); % change explicit_midpoint_step to the other funcs when needed
    error = norm(Xn1 - sol_func(t + h));
    diff = norm(sol_func(t+h) - sol_func(t));
end