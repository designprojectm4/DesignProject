function dT_rk4=dT_rk4(Fvec,Tvec,dF,dA,dW)
%remember dT=heat(Fvec,Tvec,dF,dA,dW);)
K1=heat(Fvec,Tvec,dF,dA,dW);
K2=heat(Fvec,Tvec+0.5*dW*K1,dF,dA,dW);
K3=heat(Fvec,Tvec+0.5*dW*K2,dF,dA,dW);
K4=heat(Fvec,Tvec+dW*K3,dF,dA,dW);
dT_rk4=(dW/6)*(K1+2*K2+2*K3+K4);
end