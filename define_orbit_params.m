%orbit_params: a struct describing the system parameters
% orbit_params.m_sun is the mass of the sun
% orbit_params.m_planet is the mass of the planet
% orbit_params.G is the gravitational constant
function orbit_params = define_orbit_params()
    orbit_params = struct();
    orbit_params.m_sun = 1;%1.9891e30;
    orbit_params.m_planet = 1;%5.972e24;
    orbit_params.G = 40;%9.81;

end