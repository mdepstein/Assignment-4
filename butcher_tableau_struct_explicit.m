function BT_struct = butcher_tableau_struct_explicit()
%second order
ralston = struct();
ralston.A = [0 0; 2/3 0];
ralston.B = [.25 .75];
ralston.C = [0 2/3];

%third order
ssprk3 = struct();
ssprk3.A = [0 0 0; 1 0 0;.25 .25 0];
ssprk3.B = [1/6 1/6 2/3];
ssprk3.C = [0 1 1/2];

%fourth order
RK4 = struct();
rk4.A = [0 0 0 0; .5 0 0 0; 0 .5 0 0; 0 0 1 0];
rk4.B = [1/6 1/3 1/3 1/6];
rk4.C = [0 .5 .5 1];

    BT_struct = [ralston, ssprk3, rk4];


end


