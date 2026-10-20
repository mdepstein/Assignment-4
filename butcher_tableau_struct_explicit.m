function BT_struct = butcher_tableau_struct_explicit()
%first order
forward = struct();
forward.A = [0];
forward.B = [1];
forward.C = [0];

%second order
midpoint = struct();
midpoint.A = [0 0;.5 0];
midpoint.B = [0 1];
midpoint.C = [0 .5];

%third order
kutta = struct();
kutta.A = [0 0 0; .5 0 0; -1 2 0];
kutta.B = [1/6 2/3 1/6];
kutta.C = [0 .5 1];
    BT_struct = [forward, midpoint, kutta];
end


