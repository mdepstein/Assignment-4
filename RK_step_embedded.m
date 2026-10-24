%ThisfunctioncomputesthevalueofXatthenexttimestep
%foranyarbitraryembeddedRKmethod
%INPUTS:
%rate_func_in:thefunctionusedtocomputedXdt.rate_func_inwill
% havetheform:dXdt=rate_func_in(t,X) (tisbeforeX)
%t:thevalueoftimeatthecurrentstep
%XA:thevalueofX(t)
%h:thetimeincrementforasinglestepi.e.delta_t=t_{n+1}-t_{n}
%BT_struct:astructthatcontainstheButchertableau
% BT_struct.A:matrixofa_{ij}values
% BT_struct.B:vectorofb_ivalues
% BT_struct.C:vectorofc_ivalues
%OUTPUTS:
%XB1:theapproximatevalueforX(t+h)usingthefirstrowoftheTableau
%XB2:theapproximatevalueforX(t+h)usingthesecondrowoftheTableau
%num_evals:Acountofthenumberoftimesthatyoucalled
% rate_func_inwhencomputingthenextstep
function[XB1,XB2,num_evals]=RK_step_embedded(rate_func_in,t,XA,h,BT_struct, orbit_params)
    a = BT_struct.A;
    b = BT_struct.B;
    c = BT_struct.C;
    num_evals = 1;
    k = zeros(length(b),length(XA));
    k(1,:) = rate_func_in(t,XA,orbit_params);
    for i = 2:length(b)
        j = i-1;
        k(i,:) =  rate_func_in(t + c(i)*h, XA+h*a(i,1:j)*k(1:j)', orbit_params);
        num_evals = num_evals + 1;
    end
    
    XB1 = XA + h*b(1,:)*k;
    XB2 = XA + h*b(2,:)*k;
end