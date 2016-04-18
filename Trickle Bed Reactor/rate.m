function dFvec=rate(Fvec,T,P,d_P,CSA,voidage,cat_density)
%calculates reaction rates
%Uses Kinetics from Zhou et al
%--------------------------------------------------------------------------
%C=F/Q
[~,Q_L]=VolumeFlow(Fvec,T,P);%liquid volume flow [m3/s]
%---prefactor|activation energy--------------------------------------------
k10=[01.54e4, 86.56]; %mol g-1 s-1,kJ mol-1
k20=[7.16e3, 57.80]; %mol g-1 s-1,kJ mol-1
bG0=[2.22e-3, 36.42]; %m3 mol-1,kJ mol-1
bA0=[8.73e-3, 25.94]; %m3 mol-1,kJ mol-1
bP0=[5.80e-3, 25.77]; %m3 mol-1,kJ mol-1
bH0=[1.86e-5, 36.24]; %m3 mol-1,kJ mol-1
%pre-exponential and activation energy vectors
pre=[k10(1);k20(1);bG0(1);bA0(1);bP0(1);bH0(1)];
act=[-k10(2);-k20(2);bG0(2);bA0(2);bP0(2);bH0(2)];
%Constants
Rg=8.314459848e-3;%KJ K-1 mol-1
%Selectivity to Propanols
S_POH=0.023;
%--------------------------------------------------------------------------
%Calculations
Pvec=bsxfun(@times,pre,(exp(act/(Rg*T))));
k1=Pvec(1);k2=Pvec(2);bG=Pvec(3);bA=Pvec(4);bP=Pvec(5);bH=Pvec(6);
cG=Fvec(2)/Q_L; cA=Fvec(3)/Q_L; cP=Fvec(4)/Q_L; PH=P(1);
%Remember [FH;FG;FA;FP;FW;FPr]
%Rates
rate1=k1*bG*cG/(1+bG*cG+bA*cA+bP*cP);%G=>A
rate2=k2*bA*cA*bH*PH/((1+bG*cG+bA*cA+bP*cP)*(1+sqrt(bH*PH))^2);%A+H2=>P
rate3=S_POH*rate2;%P+H2=>Propanol+H2O, percentage of rate2 assumed
%--------------------------------------------------------------------------
%For use in effectiveness factor calculation
k=rate2/cA;
    function eff = effectiveness(Fvec,T,P,d_P,CSA,voidage,cat_density,k)
    %Calculates Overall Effectiveness Factor
    %--------------------------------------------------------------------------
    %Inputs:
    fw=Burghardt(Fvec,T,P,d_P,CSA,voidage);
    %       Wetting Efficiency
    %       Particle Diameter d_P[m]
    d_P=d_P*100;%[cm]
    %       dynamic viscosity of acetol calculated using find_mu(T)
    %       mu_L=[mu_G;mu_A;mu_P;mu_W;mu_POH]
    %       rate
    %       Re_L
    [~,Re_L]=Reynold(Fvec,T,P,d_P,CSA,voidage);
    %       liquid kinematic viscosity
    [~,visc_liq]= viscosity(T,P,Fvec);
    %--------------------------------------------------------------------------
        function D_A=StokesEinstein(T)
            Tref=273+25;
            Dref=1.16e-5;%cm2/s Uses value for acetone
            [~,mu_L]=find_mu(Tref);mu_ref=mu_L(2);
            [~,mu_L]=find_mu(T);mu_A=mu_L(2);
            D_A=Dref*(T/Tref)*(mu_ref/mu_A);
        end
    D_A=StokesEinstein(T);%cm2/s
    %--------------------------------------------------------------------------
    %Liquid to Solid Mass Transfer coefficient
    %Uses correlation from Delaunay et al (1982)
    Sc=visc_liq/(D_A*1e-4);%Schmidt Number
    %uses diffusion coefficent in m2/s
    fwSh=1.84*(Re_L^0.48)*Sc^(1/3);%Sherwood number correlation
    %--------------------------------------------------------------------------
    thiele=(d_P/6)*(cat_density*k/D_A)^0.5;
    Sh_ls=fwSh/fw;
    eff=(fw/thiele)*tanh(thiele/fw)/(1+(thiele/Sh_ls)*tanh(thiele/fw));
    end
eff = effectiveness(Fvec,T,P,d_P,CSA,voidage,cat_density,k);
rate2=eff*rate2;
%--------------------------------------------------------------------------
dFG=-rate1;
dFA=rate1-rate2;
dFH=-rate2-rate3;
dFP=rate2-rate3;%I TRIED
dFW=rate3;
dFPOH=rate3;
dFvec=[dFH;dFG;dFA;dFP;dFW;dFPOH];
end