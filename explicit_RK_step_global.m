%This function computes the value of X at the next time step
%for any arbitrary RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
%   have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1}- t_{n}
%BT_struct: a struct that contains the Butcher tableau
%   BT_struct.A: matrix of a_{ij} values
%   BT_struct.B: vector of b_i values
%   BT_struct.C: vector of c_i values
% val: index indicating which Runge-Kutta method is in use
%OUTPUTS:
%XB: the approximate value for X(t+h) (the next step)
%   formula depends on the integration method used
%num_evals: A count of the number of times that you called
%   rate_func_in when computing the next step
function [XB, num_evals] = explicit_RK_step_global(rate_func_in,t,XA,h,BT_struct)
    a = BT_struct.A;
    b = BT_struct.B;
    c = BT_struct.C;
    num_evals = 1;
    XB = zeros(size(XA));
    k = zeros(length(c),length(XA));
    orbit_params = define_orbit_params();
    k(1,:) = rate_func_in(t,XA,orbit_params)';
    for i = 2:length(b)
        j = i-1;
        k(i,:) =  rate_func_in(t + c(i)*h, XA+h*a(i,1:j)*k(1:j)', orbit_params);
        num_evals = num_evals + 1;
    end
    XB = XA + h*k'*b(1,:)';


end
