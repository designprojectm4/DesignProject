function [dFA,dF_rk4] = dF_rk4(Fvec,dW,T,P,d_P,CSA,voidage,cat_density)
%RK4-Calculates Runga-Kutta Parameters-------------------------------------
%Input Parameters
%Fvec=Fvec[FH;FG;FA;FP;FW;FPr]:Component Molar Flow Rates
%dW: weight step size
%T: Reactor Temperature at dW (Might want to add another coupled equation)
%dV=dW*voidage/((cat_density)*(1-voidage)): Differential Volume Element
%--------------------------------------------------------------------------
K1=rate(Fvec,T,P,d_P,CSA,voidage,cat_density);
K2=rate(Fvec+0.5*dW*K1,T,P,d_P,CSA,voidage,cat_density);
K3=rate(Fvec+0.5*dW*K2,T,P,d_P,CSA,voidage,cat_density);
K4=rate(Fvec+dW*K3,T,P,d_P,CSA,voidage,cat_density);
dF_rk4=(dW/6)*(K1+2*K2+2*K3+K4);
dFA=dF_rk4(3);%Rate of Acetol Consumption
end
%Remember: dFvec=[dFH;dFG;dFA;dFP;dFW;dFPOH]