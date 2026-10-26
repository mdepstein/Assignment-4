%first order system driven by sinusoid
function dXdt = rate_func01(t,X,orbit_params)
    dXdt = -5*X + 5*cos(t) - sin(t);
end