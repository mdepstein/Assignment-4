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
function[XB1,XB2,num_evals]=RK_step_embedded(rate_func_in,t,XA,h,BT_struct)
    a = BT_struct.A;
    b1 = BT_struct.B(1,:);
    b2 = BT_struct.B(2,:);
    c = BT_struct.C;
    num_evals = 0;
    XB1 = zeros(size(XA));
    XB2 = zeros(size(XA));
    
    k = zeros(length(XA),length(b1));
   
    k(:,1) = rate_func_in(t,XA);

    
    for i = 2:length(b1)
        j = i-1;
        
        k(:,i) =  rate_func_in(t + c(i)*h, XA+h*k(:,1:j)*(a(i,1:j)'));
        num_evals = num_evals + 1;
    end

    XB1 = XA + h*k*b1';
    XB2 = XA + h*k*b2';
end
