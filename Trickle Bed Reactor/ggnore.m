clear
clc
nmax=1e6;
n=0;
ICS%Initial Conditions
%--------------------------------------------------------------------------
%WE ON AN ULTRALIGHT BEAM. WE ON AN ULTRALIGHT MEME. THIS IS A GOD DREAM.
%--------------------------------------------------------------------------
%Array Pre-Allocation
W_array=zeros(1e6,1);
X_array=zeros(1e6,1);
acetol_rate=zeros(1e6,1);
FH_array=zeros(1e6,1);
FG_array=zeros(1e6,1);
FA_array=zeros(1e6,1);
FP_array=zeros(1e6,1);
FW_array=zeros(1e6,1);
FPOH_array=zeros(1e6,1);
T_array=zeros(1e6,1);
dT_array=zeros(1e6,1);
Tj_array=zeros(1e6,1);
PH_array=zeros(1e6,1);
P_L_array=zeros(1e6,1);
%--------------------------------------------------------------------------
for i=1:450
n=n+1;
if n>nmax %this bit is for banter
        button=questdlg('What the fuck man?','Error','I''m an idiot','I''m stupid','I''m an idiot');
        h=msgbox('Yes you are');
        break
else
%--------------------------------------------------------------------------
%STEP CHANGES
W=W+dW;
W_array(n,:)=W;
z=z+dz;
%--------------------------------------------------------------------------
%Flow Rates
[dFA,dF] = dF_rk4(Fvec,dW,T,P);
Fvec_new=Fvec+dF;
Fvec=Fvec_new;
if Fvec(1)<=0||Fvec(3)<=0 %REMEMBER: [FH;FG;FA;FP;FW;FPOH]
    error('this is literally impossible');
end
X_A=(Fvec0(3)-Fvec(3))/Fvec0(3);
if X_A<X_Aold % u fukd up m8
    error('this is literally impossible');
end
X_Aold=X_A;
X_array(n,:)=X_A;
acetol_rate(n,:)=dFA;%rate of acetol production
FH_array(n,:)=Fvec(1);
FG_array(n,:)=Fvec(2);%
FA_array(n,:)=Fvec(3);% NOT EFFICIENT CODING
FP_array(n,:)=Fvec(4);% PRE-ALLOCATE ARRAY SIZES
FW_array(n,:)=Fvec(5);%
FPOH_array(n,:)=Fvec(6);
%REMEMBER: Fvec=[FH;FG;FA;FP;FW;FPOH]
%--------------------------------------------------------------------------
%Temperatures
dT=dT_rk4(Fvec,Tvec,dF,dA,dW);
%Calculates temperature change as a result of change in molar flows
dT_array(n,:)=dT(1);
Tvec=Tvec+dT;
T=Tvec(1);
Tj=Tvec(2);
T_array(n,:)=T;
Tj_array(n,:)=Tj;
%Pressure Drop
dP=dP_rk4(dW,Fvec,T,P,d_P,CSA,voidage,cat_density,dzdW);
P=P+dP;
PH=P(1);
PH_array(n,:)=PH;
P_L=P(2);
P_L_array(n,:)=P_L;
end
end
FG=Fvec(1);FA=Fvec(2);FP=Fvec(3);FH=Fvec(4);FW=Fvec(5);FPOH=Fvec(6);
%plot(W_array,FG_array,W_array,FA_array,W_array,FP_array,W_array,FW_array,W_array,FPOH_array)%,W_array,FH_array)
%plot(W_array,X_array)
plot(W_array,T_array);
T_max=max(T_array)
T_min=min(T_array)