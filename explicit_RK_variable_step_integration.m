%Runs numerical integration arbitrary RK method using variable time steps
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
%have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%p: how error scales with step size (error = k*hˆp)
%error_desired: the desired local truncation error at each step
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals, step_fail_rate] = explicit_RK_variable_step_integration ...
    (rate_func_in,tspan,X0,h_ref,BT_struct,p,error_desired)

    t0 = tspan(1); tf = tspan(2);
    X = X0; h = h_ref;
    t_list = zeros(1, 100);
    X_list = zeros(size(t_list));
    h_list = zeros(size(t_list));
    failure_num = 0;
    
    for i = 2:(length(t_list)-1)
        t = t_list(i);
        while t < tf
           [XB, num_evals, h_next, redo] = explicit_RK_variable_step ...
               (rate_func_in,t,X0,h,BT_struct,p,error_desired) 
           if redo == false
               X = XA;
               t_list(i+1) = t_next;
               X_list(i+1) = X;
               t = t_next;
           else
               failure_num = failure + 1;
           end
           h = h_next;
           h_list(i) = h;
        end
    end
    h_avg = mean(h_list);
    step_fail_rate = failure_num/num_evals;
end