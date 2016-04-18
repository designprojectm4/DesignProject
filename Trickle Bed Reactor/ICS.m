clear
clc
nmax=1e6;
n=0;
%--------------------------------------------------------------------------
%WE ON AN ULTRALIGHT BEAM. WE ON AN ULTRALIGHT MEME. THIS IS A GOD DREAM.
%--------------------------------------------------------------------------
%Reactor Conditions;
T=220+273; %Inlet Temperature [K]
Tj=100+273; %Inlet Cold Utility Temperature [K] Co-current
Tvec=[T;Tj];
P=3e6; %Inlet Pressure assumed for now
P=[P;P]; %[P_G(GAS);P_L(LIQUID)] 
LHSV=10.5/3600;%V=Q_L/LHSV;%VOLUME
%--------------------------------------------------------------------------
%Catalyst Bed Propertes
voidage=0.40;
cat_density=5677642.392;%g/m3
d_P=3e-3;%m
CSA=1;%Trickle Bed CSA [m2]
%uL=Q_L/CSA;%Superficial Liquid Velocity
%--------------------------------------------------------------------------
%Initial Conditions(molar flow rates in [mol/s])
FG=5.96E-28; %Glycerol
FA=8.973223392; %Acetol
FP=0;%Propylene Glycol
FH=54.53958333; %Hydrogen
%FH=300;
FW=1.087872813;%Water
%FW=18.3;%Steve's column
FPOH=0.3257064118;%Propanol
%Conversion
X_A=0;
X_Aold=0;
%Step Size
V=0;
W=0;
%FH=PH*Ng/P;
Fvec0=[FH;FG;FA;FP;FW;FPOH];%Inlet Flow Rates
%--------------------------------------------------------------------------
%Step Size
dW=10;
dV_void=dW*voidage/((cat_density)*(1-voidage));%dV is void volume
%i.e. we assume concentration=molar flow/void volume
dV=dW/((cat_density)*(1-voidage));
dz=dW/((cat_density)*(1-voidage)*CSA);
dzdW=1/((cat_density)*(1-voidage)*CSA);
z=0;
dA=2*dz*sqrt(pi*CSA);%heat transfer area
%--------------------------------------------------------------------------
Fvec=Fvec0;
%REMEMBER: [FH;FG;FA;FP;FW;FPOH]
