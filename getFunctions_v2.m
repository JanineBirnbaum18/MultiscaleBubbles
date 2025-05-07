function [SolFun, DiffFun, ViscFun, m0_fun, pb_fun,...
    PTt_fun] = getFunctions_v2(SolModel,DiffModel, ViscModel,...
    EOSModel, PTtModel)
%This script returns the functions that are used by the bubble growth
% model (Numerical_Model.m)
% 
% See the user manual for more info and references to the individual
% models.
% The references can also be found at the bottom of
% Bubble_Growth_Modelling.m
% 
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% If using or adapting this code, please cite Coumans et al., (2020):
% 
% Authors: J. P. Coumans, E. W. Llewellin*, F. B. Wadsworth
% M. C. S Humphreys, S. A.  Mathias, B. M. Yelverton, and J. E. Gardner
% 
% Title: An experimentally validated numerical model 
% for bubble growth in magma 
% 
% Journal: Journal of Volcanology and Geothermal Research (JVGR)
% 
% Year: 2020
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%Below are the switch-case statements that 'point' to the desired functions
%as defined by the user in the Bubble_Growth_Modelling.m script
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================

%==========================================================================
%Return the target Solubility function
%==========================================================================
switch SolModel
    case 'Ryan 2015'
        SolFun = @(T,P)Ryan2015_sol(T);
        
    case 'Liu 2005'
        SolFun = @(T,P)Liu2005_sol(T,P);

    case 'Weaver'
        SolFun = @(T,P)Weaver_sol(T);
       
    case 'Krafla_constant'
        SolFun = @(T,P)Krafla_sol(T);

    case 'Schunke'
        SolFun = @(T,P)Schunke_sol(T,P);
end
%==========================================================================
%Return the target Diffusivity function
%==========================================================================
switch DiffModel
    case 'Zhang 2010 Metaluminous simple'
        DiffFun = @(H2Ot,T,P,W)Zhang_2010_Metaluminous_Simple_fun(H2Ot,T,P);
        
    case 'Zhang 2010 Metaluminous'
        DiffFun = @(H2Ot,T,P,W)Zhang_2010_Metaluminous_fun(H2Ot,T,P,W);
        
    case 'Zhang 2010 Peralkaline'
        DiffFun = @(H2Ot,T,P,W)Zhang_2010_Peralkaline_fun(H2Ot,T,P,W);
        
    case 'Constant'
        DiffFun = @(H2Ot,T,P,W)ConstantDiff_fun(H2Ot);  
end

%==========================================================================
%Return the target Viscosity function
%==========================================================================
switch ViscModel
    case 'Giordano 2008'
        ViscFun = @(H2Ot,T,Composition)Giordano_2008_visc(H2Ot,T,Composition);
    
    case 'Hess and Dingwell 1996'
        ViscFun = @(H2Ot,T,Composition)HessDingwell_1996_visc(H2Ot,T);
        
    case 'Peralkaline Giordano 2000'
        ViscFun = @(H2Ot,T,Composition)PerAlkaline_Giordano_2000_visc(H2Ot,T);
end

%==========================================================================
%Return the target functions for Equation of State
%==========================================================================
switch EOSModel
    case 'Ideal Gas Law'
        m0_fun = @(R,P,T)m0_IdealGas_fun(R,P,T);
        pb_fun = @(m, T, R)pb_IdealGas_fun(m,T,R);
        
    case 'Pitzer and Sterner'
        m0_fun = @(R,P,T)m0_pitzer_fun(R,P,T);
        pb_fun = @(m, T, R)pb_Pitzer_fun(m,T,R);        
end
%==========================================================================
%Return the pressure-temperature-time functions
%==========================================================================
switch PTtModel
    case 'P: Isobaric, T: Isothermal'
        PTt_fun = @(P_0, P_f, dPdt, T_0, T_f, dTdt, t_quench, t)[P_f, T_f].*ones(size(t));
        
    case 'P: Isobaric, T: Polythermal-Dwell'
        PTt_fun = @(P_0, P_f, dPdt, T_0, T_f, dTdt, t_quench, t)...
            Isobaric_Polythermal_Dwell_fun(P_0,T_0,T_f,dTdt,t);
        
    case 'P: Polybaric-Dwell, T: Isothermal-quench'
        PTt_fun = @(P_0, P_f, dPdt, T_0, T_f, dTdt, t_quench, t)...
            Polybaric_Dwell_Isothermal_Quench_fun(P_0,P_f,dPdt,T_0,T_f,dTdt,t);  

    case 'Evolving'
        PTt_fun = @(P_0, P_f, dPdt, T_0, T_f, dTdt, t_quench, t)...
            Evolving_fun(P_0,dPdt,T_0,dTdt,t);

    case 'Quench'
        PTt_fun = @(P_0, P_f, dPdt, T_0, T_f, dTdt, t_quench, t)...
            Quench_fun(P_0,P_f,dPdt,T_0,T_f,dTdt,t_quench, t);

    case 'Experiment'
        PTt_fun = @(P_0, P_f, dPdt, T_0, T_f, dTdt, t_quench, t)...
            Experiment_fun(P_0,P_f,dPdt,T_0,T_f,dTdt,t_quench, t);

    case 'Optical dilatometer'
        PTt_fun = @(P_0, P_f, dPdt, T_0, T_f, dTdt, t_quench, t)...
            OpticalDilatometer_fun(P_0,P_f,dPdt,T_0,T_f,dTdt,t_quench,t);

    case 'Jenny'
        PTt_fun = @(P_0, P_f, dPdt, T_0, T_f, dTdt, t_quench, t)...
            Jenny_fun(P_0,P_f,dPdt,T_0,T_f,dTdt,t_quench, t);
end

%==========================================================================
%==========================================================================
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%Below are the functions that can be called by the model to perform
%Calculations. There are functions for Viscosity, Solubility, Diffusivity,
%The inital mass of gas (equation of state), and pressure in the bubble
%(equation of state)
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%==========================================================================
%==========================================================================

%==========================================================================
%Solubility functions
%==========================================================================
%Regressed solubility function from Ryan et al. (2015)
function H2Oeq = Ryan2015_sol(T)
H2Oeq = 92.3/T + 0.0287;

%Solubility function from Liu et al. (2005)
function H2Oeq = Liu2005_sol(T,P)
% convert pressure to MPa from Pa
P = P*1e-6; %Pa to MPa
H2Oeq =((354.94.*sqrt(P)+9.623.*P-1.5223.*(P.^1.5))./T) + 0.0012439.*(P.^1.5);

function H2Oeq = Weaver_sol(T)
H2Oeq = 0.088 + 0*T;

function H2Oeq = Krafla_sol(T)
H2Oeq = 0.039 + 0*T;

function H2Oeq = Schunke_sol(T,P)
P = P*1e-6;
%H2Oeq =((276.8*sqrt(P)+71.15.*P-1.5223.*(P.^1.5))./T) + 0.0012439.*(P.^1.5);
PCO2 = 0;
H2Oeq =((354.94.*sqrt(P)+9.623.*P-1.5223.*(P.^1.5))./T) + 0.0012439.*(P.^1.5) + PCO2.*(-1.084.*1E-4.*sqrt(P)-1.362.*1E-5.*P);

%==========================================================================
%Diffusivity functions
%==========================================================================
%Diffusivity function for metaluminous rhyolite from Zhang and Ni (2010)
%(Equation 15 therein)
function D = Zhang_2010_Metaluminous_Simple_fun(H2Ot,T,P)
P=P/10^9; %Convert pressure to GPa from Pa
D = H2Ot.*exp(-18.1+(1.888.*P)-((9699+3626.*P)./T));


%Diffusivity function for metaluminous rhyolite from Zhang and Ni (2010)
%(Equations 7a, 13, 14 therein)
function D = Zhang_2010_Metaluminous_fun(H2Ot,T,P,W)
P=P/10^9; %Convert pressure to GPa from Pa
%convert from wt.% to mole fraction
XH2Ot=(H2Ot/18.015)./((H2Ot/18.015)+(100-H2Ot)/W); 
%Compute the diffusivity of molecular water (Zhang & Ni 2010,
%equation 14)
DH2Om = exp(-14.26 + 1.888.*P - 37.26.*XH2Ot - ...
    ((12939 + 3626.*P - 75884.*XH2Ot)./T));
%Get the speciation constant (Zhang & Ni 2010, equation 7a)
K = exp(1.876 - 3110./T);
%Compute the total water diffusivity from Diffusivity of
%molecular water and speciation constant (Zhang & Ni 2010, equation
%13)
D = DH2Om.*(1-((0.5-XH2Ot)./((4./K-1).*(XH2Ot-(XH2Ot.^2))+0.25).^(1/2)));

%Diffusivity function for peralkaline rhyolite from Zhang and Ni (2010)
%(Equations 7a, 13, 16 therein)
function D = Zhang_2010_Peralkaline_fun(H2Ot,T,P,W)
P=P/10^9;%Convert from Pa to GPalog
%convert from wt.% to mole fraction
XH2Ot=(H2Ot/18.015)./((H2Ot/18.015)+(100-H2Ot)/W); 
%Compute the diffusivity of molecular water (Zhang & Ni 2010,
%equation 16)
DH2Om = exp(-12.79 - 27.87.*XH2Ot - ...
    ((13939 + 1230.*P - 60559.*XH2Ot)./T));
%Get the speciation constant (Zhang & Ni 2010, equation 7a)
K = exp(1.876 - 3110./T);
D = DH2Om.*(1-((0.5-XH2Ot)./((4./K-1).*(XH2Ot-(XH2Ot.^2))+0.25).^(1/2)));

%Diffusivity function for a constant diffusivty
function D = ConstantDiff_fun(H2Ot)
DiffVal = 1e-11; %Diffusivity value in m^2 s^-1
D = ones(size(H2Ot)).*DiffVal;


%==========================================================================
%Viscosity functions
%==========================================================================
%Viscosity function from Hess and Dingwell (1996)
function v = HessDingwell_1996_visc(H2Ot,T)
v=10.^(-3.545+0.833*log(H2Ot)+(9601-2368*log(H2Ot))./(T-(195.7+32.25.*log(H2Ot))));

%Viscosity function from Giordano et al. (2000)
function v = PerAlkaline_Giordano_2000_visc(H2O, T)
v=10.^(-5.9-0.286.*log(H2O)+(10775.4-394.8.*H2O)./(T-148.7+21.65.*log(H2O)));

function v = Giordano_2008_visc(H2Ot, T, Composition)
%Get the Giordano2008 VFT parameters from their script
%Giordano2008_Model.mat
outvalues=Giordano2008_Model(H2Ot, Composition);
At = outvalues(:,1);
Bt = outvalues(:,2);
Ct = outvalues(:,3);
v = 10.^(At+(Bt./(T-Ct)));


%==========================================================================
%Equation of State Initial gas mass functions
%==========================================================================
function m0 = m0_IdealGas_fun(R,P,T)
%calculate initial mass of water in the bubble (kg)
%from pV=nRT
V=(4*pi()/3)*R^3; %Volume of the bubble (M^3)
%Moles of water in the bubble (from volume and P)
%Gas constant = 8.314 J mol^-1K-1
% 1 Pa*m^3 = 1 joule
n=P*V/(8.314*T);
%18.015 g/mol * n moles = mass (g)* 1kg / 1000g
m0=18.015*n/1000;

function m0 = m0_pitzer_fun(R,P,T)
%calculate initial mass of water in the bubble (kg)
V=(4*pi()/3)*R.^3; %Volume of the bubble in M^3
%Returns the density of the bubble in kg/m^3
rho = density(P,T,coefficients());
m0 = rho.*V;

%==========================================================================
%Equation of State bubble pressure functions
%==========================================================================
%calculate internal pressure based on the ideal gas law (PV = nRT)
function pb = pb_IdealGas_fun(m, T, R)
pb=(((m.*1000)./18.015).*8.314.*T)./((4.*pi./3).*R.^3);

%calculate internal pressure based on pitzer & Sterner, 94
%Returns the pressure in the bubble for a given density (mol/cm^3),
%Temperature (kelvin) and the coefficient matrix
function pb = pb_Pitzer_fun(m, T, R)

%Volume of the bubble in m^3
Vbub = (4/3)*pi().*R .^3;

%Rho must be in terms of mol/cm^3 for the pitzer & sterner equation of
%state. Therefore Volume: m^3 * (1e6 cm^3 / 1 m ^3.
%and Mass: (kg * (1000g / 1 kg))*(1 mole / 18.015 g)
rho = ((m.*1000)./18.015)./(Vbub.*1e6);

%Get the pitzer and sterner coefficients
b = coefficients();
col = size(T,2);

a=zeros(10,col);
for i=1:10
    a(i,:)=b(i,1).*T.^-4 + b(i,2).*T.^-2 + b(i,3).*T.^-1 +...
        b(i,4) + b(i,5).*T + b(i,6).*T.^2;
end

%P [bars], R [cm^3*bar/K/mol] and T [K]

%NOTE Gas constant is in terms of [cm^3*bar/K/mol] as the equation of PS,
%94 has pressure in terms of bar and density in terms of mol cm^-3
pb = (rho+a(1,:).*rho.^2-rho.^2.*((a(3,:)+2.*a(4,:).*rho+3.*a(5,:).*rho.^2+4.*a(6,:).*rho.^3)./((a(2,:)+a(3,:).*rho+a(4,:).*rho.^2+a(5,:).*rho.^3+a(6,:).*rho.^4).^2))+a(7,:).*rho.^2.*exp(-a(8,:).*rho)+a(9,:).*rho.^2.*exp(-a(10,:).*rho)).*(83.14472.*T);

% Convert P from bars (equation P&S) to pascals (model)
pb = pb/1e-5; %1 pascal a 1e-5 bars



%==========================================================================
%P-T-t functions
%==========================================================================
function out = Isobaric_Polythermal_Dwell_fun(P_0,T_0,T_f,dTdt,t)
T_t_ramp = abs((T_f-T_0)/dTdt);
T = ((t<=T_t_ramp).*(T_0 + t.*dTdt))+((t>T_t_ramp).*T_f);
P = P_0.*ones(size(t));
out = [P, T];

function out = Polybaric_Dwell_Isothermal_Quench_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t)
T_t_Quench = abs((T_f-T_0)/dTdt);
P_t_Press = abs((P_f-P_0)/dPdt);
P = ((t<=P_t_Press).*(P_0 + t.*dPdt))+((t>P_t_Press).*P_f);
T = ((t<=P_t_Press).*T_0) + ...
    (((t>P_t_Press) & (t<=P_t_Press+T_t_Quench)).*(T_0 + (t-P_t_Press).*dTdt)) + ...
    +((t>P_t_Press+T_t_Quench).*T_f);
out = [P, T];

function out = Evolving_fun(P_0,dPdt,T_0,dTdt,t)
P = P_0 + t.*dPdt;% + t.*dPdt;
T = T_0 + t.*dTdt;% + t.*dTdt;
out = [P,T];

function out = Quench_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t_quench,t)
T_t_Quench = abs((T_f-T_0)/dTdt);
P_t_Press = abs((P_f-P_0)/dPdt);
P = ((t<=P_t_Press).*(P_0 + t.*dPdt))+((t>P_t_Press).*P_f);
T = ((t<=t_quench).*T_0) + ...
    (((t>t_quench) & (t<=t_quench+T_t_Quench)).*(T_0 + (t-t_quench).*dTdt)) + ...
    +((t>t_quench+T_t_Quench).*T_f);
out = [P, T];

function out = Experiment_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t_quench,t)
T_t_hold1 = abs((T_f(1)-T_0)/dTdt);
T_t_hold2 = T_t_hold1 + t_quench(1) + abs((T_f(2)-T_f(1))/dTdt);
P_t_hold = abs((P_f-P_0)/dPdt);
P = ((t<=P_t_hold).*(P_0 + t.*dPdt))+((t>P_t_hold).*P_f);
T = (t<=T_t_hold1).*(T_0 + t.*dTdt) + ...
    ((t>T_t_hold1) & (t<=T_t_hold1 + t_quench(1))).*(T_f(1)) + ...
    ((t>T_t_hold1 + t_quench(1)) & (t<=T_t_hold2)).*(T_f(1) + (t-(T_t_hold1+t_quench(1))).*dTdt) + ...
    ((t>T_t_hold2) & (t<=T_t_hold2+t_quench(2))).*(T_f(2)) + ...
    (t>T_t_hold2 + t_quench(2)).*(T_f(2) - 20*sqrt(abs(t-(T_t_hold2 + t_quench(2)))/30+5) + 20*sqrt(5));
out = [P, T];

function out = OpticalDilatometer_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t_quench,t)
T_t_hold1 = abs((T_f(1)-T_0)/dTdt);
T_t_hold2 = T_t_hold1 + t_quench(1) + abs((T_f(2)-T_f(1))/dTdt);
P_t_hold = abs((P_f-P_0)/dPdt);
P = ((t<=P_t_hold).*(P_0 + t.*dPdt))+((t>P_t_hold).*P_f);
T = (t<=T_t_hold1).*(T_0 + t.*dTdt) + ...
    ((t>T_t_hold1) & (t<=T_t_hold1 + t_quench(1))).*(T_f(1)) + ...
    ((t>T_t_hold1 + t_quench(1)) & (t<=T_t_hold2)).*(T_f(1) + (t-(T_t_hold1+t_quench(1))).*dTdt) + ...
    ((t>T_t_hold2) & (t<=T_t_hold2+t_quench(2))).*(T_f(2)) + ...
    (t>T_t_hold2 + t_quench(2)).*(T_f(2) - 2/60*(t-T_t_hold2-t_quench(2)));
out = [P, T];

function out = Jenny_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t_quench,t)
T_t_Quench = abs((T_f-T_0)/dTdt(1));
P_t_Press = abs((P_f-P_0)/dPdt);
P = ((t<=P_t_Press).*(P_0 + t.*dPdt))+((t>P_t_Press).*P_f);
T = (t<=T_t_Quench).*(T_0 + dTdt(1)*t) + ...
    ((t>T_t_Quench) & (t<=t_quench+T_t_Quench)).*(T_f) + ...
    +((t>t_quench+T_t_Quench).*(T_f+dTdt(2)*(t-t_quench-T_t_Quench)));
out = [P, T];


%==========================================================================
%==========================================================================
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%Below are the additional functions that are used by the various model
%functions for performing calculations
%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
%==========================================================================
%==========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The below functions are used for solving equation of
%state via Pitzer & Sterner, 1994
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The coefficients from Pitzer and Sterner, 1994
function Coeff = coefficients()
% matrix of coefficients for eqs. in Pitzer & Sterner 1994
b=zeros(10,6);
b(1,3)=0.24657688e6;
b(1,4)=0.51359951e2;
b(2,3)=0.58638965e0;
b(2,4)=-0.28646939e-2;
b(2,5)=0.31375577e-4;
b(3,3)=-0.62783840e1;
b(3,4)=0.14791599e-1;
b(3,5)=0.35779579e-3;
b(3,6)=0.15432925e-7;
b(4,4)=-0.42719875e0;
b(4,5)=-0.16325155e-4;
b(5,3)=0.56654978e4;
b(5,4)=-0.16580167e2;
b(5,5)=0.76560762e-1;
b(6,4)=0.10917883e0;
b(7,1)=0.38878656e13;
b(7,2)=-0.13494878e9;
b(7,3)=0.30916564e6;
b(7,4)=0.75591105e1;
b(8,3)=-0.65537898e5;
b(8,4)=0.18810675e3;
b(9,1)=-0.14182435e14;
b(9,2)=0.18165390e9;
b(9,3)=-0.19769068e6;
b(9,4)=-0.23530318e2;
b(10,3)=0.92093375e5;
b(10,4)=0.12246777e3;
Coeff = b;

%Pitzer and Sterner, 1994 gas density function, calls PS_myfun
function rho = density(P,T,b)
% convert P to bars from input (pascals)
P = P.*1e-5; %1 pascal a 1e-5 bars

a=zeros(10);
for i=1:10
    a(i)=b(i,1)*T^-4 + b(i,2)*T^-2 + b(i,3)*T^-1 +...
        b(i,4) + b(i,5)*T + b(i,6)*T^2;
end
% PRT = P/RT where P [bars], R [cm^3*bar/K/mol] and T [K]
PRT = P/(83.14472*T);

% solve implicit equation for rho and convert to kg/m^3
rho = fzero(@PS_myfun,0.001,[],a,PRT)*18.01528*1000;

% the function from Pitzer & Sterner 1994, which takes the matrix of
% coefficients a and P/RT as arguments; rho is a first guess for the
% density [g/mol]
function y = PS_myfun(rho,a,PRT)
y = (rho+a(1)*rho^2-rho^2*((a(3)+2*a(4)*rho+3*a(5)*rho^2+4*a(6)*rho^3)/((a(2)+a(3)*rho+a(4)*rho^2+a(5)*rho^3+a(6)*rho^4)^2))+a(7)*rho^2*exp(-a(8)*rho)+a(9)*rho^2*exp(-a(10)*rho)) - PRT;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The below functions are used for solving viscosity
%using Giordano 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outvalues=Giordano2008_Model(H2Ot, Composition)
%This is code is from Girodano 2008 (see citation below) modified
%for use by the Coumans et al., 2020 bubble growth model

% SCRIPT grdmodel08 Nov 2007
%
% MATLAB script to compute silicate melt viscosity.
%
% Citation: Giordano D, Russell JK, & Dingwell DB (2008) 
%  Viscosity of Magmatic Liquids: A Model. Earth & Planetary Science
%  Letters, v. 271, 123-134.
%
% ________________________________________________________________
%INPUT: Chemical compositions of silicate melts (Wt. % oxides) as:
% SiO2 TiO2 Al2O3 FeO(T) MnO MgO CaO Na2O K2O P2O5 H2O F2O-1
% One line for each melt composition
% _________________________________________________________________ 

% ________________________________________________________________
% OUTPUT: VFT Model parameters A, B and C for each melt composition
%         to model temperature dependence: log n = A + B/[T(K) - C]
% VFT_model.out contains: No. , A_value , B_value , C_value , Tg , F
% VFT_curves.out contains: Temperature Array in Kelvins &
%                          Log(n) values at T(K) values (1 line per melt) 
% ________________________________________________________________

% VFT Multicomponent-Model Coefficients
% -----------------------------------------------------------------
AT      = -4.55;
bb  = [159.56  -173.34 72.13 75.69 -38.98 -84.08 141.54 -2.43 -0.91 17.62];
cc  = [2.75 15.72 8.32 10.2 -12.29 -99.54 0.3 ];

row = length(H2Ot);
Comp_rep = repmat(Composition,row,1);
Comp_rep(:,11) = H2Ot';

wtm=Comp_rep;

[nx,nc]=size(wtm);

bb_rep = repmat(bb,nx,1);
cc_rep = repmat(cc,nx,1);
AT_rep = repmat(AT,nx,1);

% Function molefrac_grd: converts wt % oxide bais to mole % oxide basis
%[ncomps xmf_t] = molefrac_grd(wtm);
[ncomps, xmf_t] = molepct_grd(wtm);
             
%outvalues=[];
% Load composition-basis matrix for multiplication against model-coefficients
% Result is two matrices bcf[nx by 10] and ccf[nx by 7]
    siti    =   xmf_t(:,1) + xmf_t(:,2);
    tial    =   xmf_t(:,2)+xmf_t(:,3);
    fmm     =   xmf_t(:,4) + xmf_t(:,5) + xmf_t(:,6);
    nak     =   xmf_t(:,8) + xmf_t(:,9);
    b1  =   siti;
    b2  =   xmf_t(:,3);
    b3  =   xmf_t(:,4) + xmf_t(:,5) + xmf_t(:,10);
    b4  =   xmf_t(:,6);
    b5  =   xmf_t(:,7);
    b6  =   xmf_t(:,8) + xmf_t(:,11) + xmf_t(:,12);
    b7  =   xmf_t(:,11) + xmf_t(:,12) + log(1+xmf_t(:,11));
    b12 =   siti.*fmm;
    b13 =   (siti + xmf_t(:,3) + xmf_t(:,10)).*( nak + xmf_t(:,11) );
    b14 =   xmf_t(:,3).*nak;

    c1      =   xmf_t(:,1);
    c2      =   tial; 
    c3      =   fmm;
    c4      =   xmf_t(:,7);
    c5      =   nak;
    c6      =   log(1+xmf_t(:,11) + xmf_t(:,12));
    c11     =   xmf_t(:,3) + fmm + xmf_t(:,7) - xmf_t(:,10);
    c11     =   c11.*(nak + xmf_t(:,11) + xmf_t(:,12)); 
    bcf      =   [b1 b2 b3 b4 b5 b6 b7 b12 b13 b14];
    ccf      =   [c1 c2 c3 c4 c5 c6 c11];   
    
%     for iz = 1:nx                     % step through each composition 
%     BT          = sum(bb.*bcf(iz,:));
%     CT          = sum(cc.*ccf(iz,:));
%     TG          = BT/(12-AT) + CT;
%     F           = BT/(TG*(1 - CT/TG)*(1 - CT/TG));
%     outvalues   =[outvalues ; iz AT BT CT TG F];
%     end
%     
%Calculate the coefficients using matrix algebra instead of a loop
BT          = sum(bb_rep.*bcf(:,:),2);
CT          = sum(cc_rep.*ccf(:,:),2);
TG          = BT./(12-AT_rep) + CT;
F           = BT./(TG.*(1 - CT./TG).*(1 - CT./TG));

outvalues   =[AT_rep BT CT TG F];

%Function for computing the mole percent
function [nox xmf] = molepct_grd(wtm)
[nr nc] =   size(wtm);
% [n x]=MOLEPCT(X) Calculates mole percent oxides from wt % 
% 1 SiO2 2 TiO2  3 Al2O3  4 FeO  5 MnO  6 MgO  7CaO 8 Na2O  9 K2O  10 P2O5 11 H2O 12 F2O-1  
% Output: mole fractions of equivalent
mw=[60.0843, 79.8658, 101.961276, 71.8444, 70.937449,40.3044,56.0774, 61.97894, 94.1960, 141.9446,18.01528, 18.9984];
mp=[];
xmf=[];

%Replicate the molecular weight row by the number of spatial nodes used
%(number of data rows)
mw_rep = repmat(mw,nr,1);

%Perform the calculation using matrix algebra
wtn = [wtm(:,1:10).*(100-wtm(:,11))./(sum(wtm(:,1:10),2)+wtm(:,12)) wtm(:,11) 0.5.*wtm(:,12).*(100-wtm(:,11))./(sum(wtm(:,1:10),2)+wtm(:,12))];
mp=wtn./mw_rep;
mpv= 100*(mp./sum(mp,2));
xmf=mpv;

[nr, nc]=size(xmf);
nox =nc;
    


%==========================================================================
% References
%==========================================================================
% Giordano, D., Dingwell, D.B. and Romano, C., 2000. Viscosity of a 
% Teide phonolite in the welding interval. J. Volcanol. Geotherm. Res., 103(1): 239-245.
% 
% Giordano, D., Russell, J.K. and Dingwell, D.B., 2008. Viscosity of
% magmatic liquids: A model. Earth Planet. Sci. Lett., 271(1): 123-134.
% 
% Hess, K.-U. and Dingwell, D.B., 1996. Viscosities of hydrous
% leucogranitic melts: A non-Arrhenian model. Am. Mineral., 81(9-10): 1297-1300.
% 
% Iacono-Marziano, G., Schmidt, B.C. and Dolfi, D., 2007. Equilibrium and
% disequilibrium degassing of a phonolitic melt (Vesuvius AD 79 “white pumice”) 
% simulated by decompression experiments. J. Volcanol. Geotherm. Res., 161(3): 151-164.
% 
% Liu, Y., Zhang, Y. and Behrens, H., 2005. Solubility of H2O in rhyolitic
% melts at low pressures and a new empirical model for mixed H2O–CO2 solubility 
% in rhyolitic melts. J. Volcanol. Geotherm. Res., 143(1): 219-235.
% 
% Mangan, M. and Sisson, T., 2000. Delayed, disequilibrium degassing 
% in rhyolite magma: decompression experiments and implications for
% explosive volcanism. Earth Planet. Sci. Lett., 183(3-4): 441-455.
% 
% Mourtada-Bonnefoi, C.C. and Laporte, D., 2004. Kinetics of bubble
% nucleation in a rhyolitic melt: an experimental study of the effect of
% ascent rate. Earth Planet. Sci. Lett., 218(3-4): 521-537.
% 
% Pitzer, K.S. and Sterner, S.M., 1994. Equations of state valid
% continuously from zero to extreme pressures for
% H2O and CO2. J. Chem. Phys., 101(4): 3111-3116.
% 
% Ryan, A.G., Russell, J.K., Nichols, A.R., Hess, K.-U. and
% Porritt, L.A., 2015. Experiments and models on H2O retrograde
% solubility in volcanic systems. Am. Mineral., 100(4): 774-786.
% 
% Tuffen, H. and Castro, J.M., 2009. The emplacement of an obsidian
% dyke through thin ice: Hrafntinnuhryggur, Krafla Iceland.
% J. Volcanol. Geotherm. Res., 185(4): 352-366.
% 
% Zhang, Y. and Ni, H., 2010. Diffusion of H, C, and O components
% in silicate melts. Rev. Mineral. Geochem., 72(1): 171-225.
































