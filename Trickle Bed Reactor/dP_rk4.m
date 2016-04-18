function dP_rk4 = dP_rk4(dW,Fvec,T,P,d_P,CSA,voidage,cat_density,dzdW)
%RK4-Calculates Runga-Kutta Parameters for Pressure Drop
K1=Holub(Fvec,T,P,d_P,CSA,voidage,cat_density,dzdW);
K2=Holub(Fvec,T,P+0.5*dW*K1,d_P,CSA,voidage,cat_density,dzdW);
K3=Holub(Fvec,T,P+0.5*dW*K2,d_P,CSA,voidage,cat_density,dzdW);
K4=Holub(Fvec,T,P+dW*K3,d_P,CSA,voidage,cat_density,dzdW);
dP_rk4=(dW/6)*(K1+2*K2+2*K3+K4);
end
%remember dP=[dP_G;dP_L]