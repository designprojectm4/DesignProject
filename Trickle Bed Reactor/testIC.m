clear
clc
T=493;
P=3e6;
FG=5.96E-28; %Glycerol
FA=8.973223392; %Acetol
FP=0;%Propylene Glycol
FH=54.53958333; %Hydrogen
FW=1.087872813;%Water
FPOH=0.3257064118;%Propanol
Fvec=[FH;FG;FA;FP;FW;FPOH];
clearvars -except T P Fvec