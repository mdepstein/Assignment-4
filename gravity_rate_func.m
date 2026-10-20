%Rate function describing Newton's law of gravitation
%INPUTS:
%t: the current time
%V: the vector of the position and velocity of the planet
% V = [x_p; y_p; dxdt_p; dydt_p]
%orbit_params: a struct describing the system parameters
% orbit_params.m_sun is the mass of the sun
% orbit_params.m_planet is the mass of the planet
% orbit_params.G is the gravitational constant
%OUTPUTS:
%dVdt: a column vector describing the time derivative of V:
% dVdt = [dxdt_p; dydt_p; d2xdt2_p; d2ydt2_p]
function dVdt = gravity_rate_func(t,V,orbit_params)
x = V(1);
y = V(2);
dxdt_p = V(3);
dydt_p = V(4);
G = orbit_params.G;
mp = orbit_params.m_planet;
ms = orbit_params.m_sun;

d2xdt2_p = -G*ms*x/(x^2+y^2)^(3/2);
d2ydt2_p = -G*ms*y/(x^2+y^2)^(3/2);
dVdt = [dxdt_p; dydt_p; d2xdt2_p; d2ydt2_p];
end