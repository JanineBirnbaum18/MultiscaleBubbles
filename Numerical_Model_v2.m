function [t, R, phi, P, T, x_out, H2Ot_all, Nb, pb_out, m] =  Numerical_Model_v2(Composition,...
    SolModel, DiffModel, ViscModel, EOSModel,OutgasModel,...
    PTtModel, SurfTens, melt_Rho, Nodes, R_0, H2Ot_0, m_loss, Nb, t_nuc, t_f, ...
    T_0, T_f, dTdt, P_0, P_f, dPdt, t_quench, Numerical_Tolerance, eta, z_p, j, Geometry, radius)
%This script defines the algorithm that solves the bubble-growth
%numerical experiment outlined in Bubble_Growth_Modelling.m
% 
% 
% See the user manual for more info.
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

%Get the user selected functions used in the numerical model
[SolFun, DiffFun, ViscFun, m0_fun, pb_fun,...
    PTt_fun] = getFunctions_v2(SolModel,DiffModel, ViscModel,...
    EOSModel, PTtModel);

%Get the molar mass of anhydrous  melt on a single oxygen basis
W = Mass_SingleOxygen(Composition);

%Get the initial porisity of the shell model
Nb_0 = Nb;
phi_0 = Porosity (Nb,R_0);

%distance to shell edge (m) (computed from R a v n
  % d phi)
L = R_0*(phi_0^(-1/3) - 1);

%Set the numerical duration bounds (from nucleation to end)
tspan= [t_nuc t_f];

%Get the P and T at the defined t_0 using the P-T-t pathway
PT_t0 = PTt_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t_quench,t_nuc);
P_t0 = PT_t0(:,1);
T_t0 = PT_t0(:,2);

%generate logarithmically spaced xB values
%scale to appropriate length, and shift so that xB
%starts at the bubble wall
%(equation A.6 from manuscript)
xB=[0 logspace(-2,0,Nodes)*L]' + R_0;


%block-centered grid - use midpoints of xB
%(equation A.3 from manuscript)
x=(xB(2:end,1)+xB(1:end-1,1))/2;

%uniform concentration; system in equilibrium Although this can be
%changed to represent non-uniform conditions. Any initial concentration
%profile can be set by the user.
if length(H2Ot_0) == 1
    H2Ot = zeros(Nodes,1)+H2Ot_0;
else
    H2Ot = H2Ot_0(:,2);
    H2Ot_0 = H2Ot_0(:,2);
end

%Compute the initial mass of gas in the bubble
if length(H2Ot_0) == 1
    m_0 = m0_fun(R_0, P_t0 +(2*SurfTens/R_0), T_t0);
else
    m_0 = m_loss;
end

%Y_0 is a vector containing quantities to be solved for by ode15s
Y_0 = [H2Ot;R_0];

%Declare the listener function
xoverFcn=@(t,X)eventDetection(t,X);

%Set the options for the ODE solver
options = odeset('AbsTol',Numerical_Tolerance(1),'RelTol',Numerical_Tolerance(2),'Events',xoverFcn);

%Create an anonymous function of the ODE to solve for passing to the
%solver
funsolve = @(t,Y)MYodeFun(t,Y,x,xB,m_0,melt_Rho,T_0,P_0,H2Ot_0,R_0,W,SurfTens,...
    SolFun,DiffFun,ViscFun,pb_fun,OutgasModel,PTt_fun,Composition,Nb,t_f,...
    T_f, dTdt, P_f, dPdt,t_quench, eta, z_p, j, Geometry,radius);

tic;

%Solve the PDE by method of lines using the ODE solver function ODE15s
[t,Y]=ode15s(funsolve,[t_nuc t_f],Y_0,options);

pb_out = 0*t;
for i = 1:length(t)
    [~,pb_out(i)] = funsolve(t(i),Y(i,:)');
end

%Get the outputs
[R, phi, P, T, x_out, H2Ot_all, Nb, m] = Outputs(Nodes,R_0,L,Y,t,Nb_0,Nb,H2Ot_0,melt_Rho,PTt_fun,P_0, P_f, dPdt,T_0,T_f,dTdt,t_quench,m_0,x);

%==========================================================================
%ODE function to be solved. See appendix A for an explanation and
%pseudocode
%==========================================================================
function [dYdt, pb] =  MYodeFun(t,Y,x,xB,m_0,melt_Rho,T_0,P_0,H2Ot_0,R_0,W,SurfTens,...
    SolFun,DiffFun,ViscFun,pb_fun,OutgasModel,PTt_fun,Composition,Nb,t_f,...
    T_f, dTdt, P_f, dPdt,t_quench, eta, z_p, j, Geometry,radius)

%extract individual concentrations
nx = (size(Y,1) - 1);
H2Ot = Y(1:nx,:);
R = Y(end,:);

%Get current temperature and pressure
PT = PTt_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t_quench,t);
P = PT(:,1);
T = PT(:,2);

%Trapz is a function that uses the trapezoidal rule to numerically
%integrate. Dim is the dimension to operate on
dim = 1;

%Get mass of gas (Equation 7 from the main manuscript)
I1=trapz(x.^3,(1/100)*H2Ot_0,dim);
I2=trapz(x.^3,(1/100)*H2Ot,dim);                                                                                                                                                                                                                           
m=m_0+4.*pi.*melt_Rho.*(I1-I2);

%Compute pressure of the gas  in the bubble from EOS
%(section 2.3.2 of main manuscript)
pb = max([1e5,pb_fun(m, T, R)]);

%====boundary conditions====
%Determine the solubility condition of water in the system based on gas
%pressure in the bubble
%(section 2.3.1 of main manuscript)
H2Oeq = SolFun(T,pb);

%get diffusion coefficients at boundaries
%Creates a vector where the first value is the boundary condition as
%determined from the equilibrium solubility of the system.
%(section 2.3.1 of main manuscript)
H2O_BC= [H2Oeq ; H2Ot];
DH2Ot = DiffFun(H2O_BC,T,P,W);

%====Solve water diffusion==== (equation A.4 from manuscript)
%molecular diffusion flux

z = [xB(1); x]; 

h1 = [z(2)-z(1); z(2:end-1)-z(1:end-2); z(end-1)-z(end-2)];
h2 = [z(3)-z(2); z(3:end)-z(2:end-1); z(end)-z(end-1)];
A = (2*h1 + h2)./h1./(h1+h2);
B = (h1+h2)./h1./h2;
C = h1./(h1+h2)./h2;
D = h2./h1./(h1+h2);
E = (h1-h2)./h1./h2; 
F = (h1 + 2*h2)./h2./(h1+h2);

dH2Otdx = [-A(1)*H2O_BC(1) + B(1)*H2O_BC(2) - C(1)*H2O_BC(3); ...
        -D(2:end-1).*H2O_BC(1:end-2) - E(2:end-1).*H2O_BC(2:end-1) + C(2:end-1).*H2O_BC(3:end); ...
        D(end).*H2O_BC(end-2) - B(end).*H2O_BC(end-1) + F(end).*H2O_BC(end)];
d2H2Otdx2 = (2*h2(2:end-1)./(h1(2:end-1).*h2(2:end-1).*(h1(2:end-1)+h2(2:end-1)))).*H2O_BC(1:end-2) +...
    -2*(h1(2:end-1)+h2(2:end-1))./(h1(2:end-1).*h2(2:end-1).*(h1(2:end-1)+h2(2:end-1))).*H2O_BC(2:end-1) + ...
    2*h1(2:end-1)./(h1(2:end-1).*h2(2:end-1).*(h1(2:end-1)+h2(2:end-1))).*H2O_BC(3:end);
d2H2Otdx2 = [d2H2Otdx2(1); d2H2Otdx2; d2H2Otdx2(end)];

dDH2Otdx  = [-A(1)*DH2Ot(1) + B(1)*DH2Ot(2) - C(1)*DH2Ot(3); ...
        -D(2:end-1).*DH2Ot(1:end-2) - E(2:end-1).*DH2Ot(2:end-1) + C(2:end-1).*DH2Ot(3:end); ...
        D(end).*DH2Ot(end-2) - B(end).*DH2Ot(end-1) + F(end).*DH2Ot(end)];

dJH2Odx = DH2Ot.*d2H2Otdx2 + dDH2Otdx.*dH2Otdx + 2./z.*DH2Ot.*dH2Otdx;

%JH2Ot = -DH2Ot.*diff([H2Oeq;H2Ot],1,1)./diff([xB(1,:);x],1,1);
%Gradient of the diffusion flux.
%dJH2Odx = (1./(x.^2)).*diff([(((x.^3+R.^3-R_0.^3).^(4/3))./(x.^2)).*JH2Ot;0],1,1)./diff(xB,1,1);

%====solve hydrodynamic equation====
%Compute the viscosity
v = ViscFun(H2Ot,T,Composition);

%Compute integrated viscosity (equation A.5 from manuscript)
I3=trapz(x, (v.*x.^2)./((R.^3-R_0.^3+x.^3).^2),dim);
%switch Geometry
%     case 'Radial'
     %I3_susp = trapz(z_p, (eta.*z_p.^2)./((z_p.^3).^2))./(4/3.*pi*(z_p(2).^3-z_p(1).^3).*Nb);
%     z_p_all = [0, z_p, radius]; 
%     dz_p = (4/3.*pi*(z_p_all(j+1).^3-z_p_all(j).^3))*(R.^3-R_0.^3)*Nb;
%     I3_right = 0;
%     if j < length(z_p)
%         z_p0 = [z_p((j+1):end), radius] - z_p(j);
%         I3_right = trapz(z_p0, ([eta((j+1):end) eta(end)].*z_p0.^2)./(((z_p(j) + dz_p).^3-z_p(j).^3 + z_p0.^3).^2))./(4/3.*pi*(z_p_all(j+1).^3-z_p_all(j).^3).*Nb);
%     end
 
%     dz_p = (4/3.*pi*(z_p_all(j+2).^3-z_p_all(j+1).^3))*(R.^3-R_0.^3)*Nb;
%     I3_left = 0;
%     if j > 1
%         z_p0 = [0, z_p(1:j-1)] - z_p(j);
%         I3_left = trapz(z_p0, ([eta(1) eta(1:j-1)].*z_p0.^2)./(((z_p(j) + dz_p).^3-z_p(j).^3 + z_p0.^3).^2))./(4/3.*pi*(z_p_all(j+2).^3-z_p_all(j+1).^3).*Nb);
%     end
     
     %I3 = I3 + min([I3_right,I3_left]);
 
 %   case 'Cylindrical'
 %   if length(z_p)>1
 %        z_p0 = [z_p, radius] - z_p(1);
 %        I3=I3 + 1./(12.*R.^2)./radius.*trapz(z_p0(2:end), ([eta(2:end) eta(end)]./z_p0(2:end).^4))./(pi.*radius.^2.*(z_p0(2) - z_p0(1)).*Nb);
 %   else
 %        I3=I3 + 1./(12.*R.^2)./radius*1/2*(eta./radius.^5)./(pi.*radius.^2.*(radius - z_p).*Nb);% + 1./(12.*R.^2)./radius.*(eta./z_p.^4)*(z_p/2)./(pi.*radius.^2.*(z_p).*Nb);
 %    end
%end

%Solve Rayleigh-Plesset equation (equation A.6 from manuscript)
dRdt= ((pb-P-(2.*(SurfTens)./R))./(12.*R.^2))./I3;

%return rhs of ode
dYdt = real([dJH2Odx(2:end);dRdt]);
switch OutgasModel
    case 'Diffusive'
        dYdt(end-1) = 0;
end

%figure(5); hold on;
%plot(x,H2Ot)
%drawnow


%==========================================================================
%Functions to get outputs from the ODE solution
%==========================================================================
function [R, phi, P, T, x_out, H2Ot_all, Nb, m] = Outputs(Nodes,R_0,L,Y,t,Nb_0,Nb,H2Ot_0,melt_Rho,PTt_fun,P_0, P_f, dPdt,T_0,T_f,dTdt,t_quench,m_0,x)
Y = Y';
%t = t';

%Get the bubble radius and phi
R = Y(end,:);
%Get all of the water profiles
H2Ot_all = Y(1:end-1,:);
%phi = Porosity (Nb,R);
[phi,Nb] = Porosity_conc (Nb_0,R_0,R);

%Get the shell thickness for each R(t)
x_out = (x.^3+R.^3-R_0^3).^(1/3);

I1=trapz(x_out(:,1),(1/100)*H2Ot_all(:,1).*x_out(:,1).^2);
I2=trapz(x_out(:,end),(1/100)*H2Ot_all(:,end).*x_out(:,end).^2);
m=m_0+4.*pi.*melt_Rho.*(I1-I2);

%Get P-T-t history
PT = PTt_fun(P_0, P_f, dPdt,T_0,T_f,dTdt,t_quench,t);
P = PT(:,1);
T = PT(:,2);

%==========================================================================
%Functions that return constants which are used in the modelling
%==========================================================================

%This function returns the molar mass of the anhydrous melt on a single
%oxygen basis
function W = Mass_SingleOxygen(Composition)

comp = Composition';

%Convert composition matrix from Viscosity input to Shishkina format
%!! SiO2 TiO2 Al2O3 FeO(T) MnO MgO CaO Na2O K2O P2O5 H2O F2O-1 !!
%Convert To
%!! SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 Cr2O3] !!

X = zeros(12,1);
X(1) = comp(1);
X(2) = comp(2);
X(3) = comp(3);
X(4) = 0;
X(5) = comp(4);
X(6) = comp(5);
X(7) = comp(6);
X(8) = comp(7);
X(9) = comp(8);
X(10) = comp(9);
X(11) = comp(10);
X(12) = 0;


Total_Mass = sum(X);

%Molar mass (g/mol) of individual elements
mSi = 28.0855;
mTi = 47.867;
mAl = 26.981539;
mFe = 55.845;
mMn = 54.938044;
mMg = 24.305;
mCa = 40.078;
mNa = 22.989769;
mK = 39.0983;
mP = 30.973762;
mCr = 51.9961;
mO = 15.999;

% [SiO2 TiO2 Al2O3 Fe2O3 FeO MnO MgO CaO Na2O K2O P2O5 Cr2O3]
%Molar mass (g/mol) of oxides
OxideMolarMass = zeros(12,1);
OxideMolarMass(1) =  (mSi+2*mO);
OxideMolarMass(2) = (mTi+2*mO);
OxideMolarMass(3) = (2*mAl+3*mO);
OxideMolarMass(4) = (2*mFe+3*mO);
OxideMolarMass(5) = (1*mFe+1*mO);
OxideMolarMass(6) = (1*mMn+1*mO);
OxideMolarMass(7) = (1*mMg+1*mO);
OxideMolarMass(8) = (1*mCa+1*mO);
OxideMolarMass(9) = (2*mNa+1*mO);
OxideMolarMass(10) = (2*mK+1*mO);
OxideMolarMass(11) = (2*mP+5*mO);
OxideMolarMass(12) = (2*mCr+3*mO);

%Compute number of moles of element, and Cation Fraction
numMolesOxygen = [2 2 3 3 1 1 1 1 1 1 5 3]';
numMolesElement = [1 1 2 2 1 1 1 1 2 2 2 2]';

%Compute the number of moles of each oxide
Moles_Oxide = X./OxideMolarMass;

%Compute moles of oxygen by stoichiometry
Moles_Oxygen = Moles_Oxide.*numMolesOxygen;

%W_melt is the mass of anhydrous melt per mole of oxygen
W = Total_Mass./sum(Moles_Oxygen);

%Compute the gas volume fraction from Nb and radius
function phi = Porosity (Nb,R)
meltvolume = 1/Nb-(4/3).*pi()*R.^3;
gasvolume = (4/3).*pi()*R.^3;
phi = gasvolume./(gasvolume+meltvolume);

function [phi,Nb] = Porosity_conc (Nb_0,R_0,R)
meltvolume = 1/Nb_0-(4/3).*pi()*R_0.^3;
gasvolume = (4/3).*pi()*R.^3;
phi = gasvolume./(meltvolume + gasvolume);
%Nb = 1./((4/3).*pi()*R.^3 + 1./Nb_0);
%Nb = Nb_0*(1 + 4/3*pi*R_0.^3.*Nb_0)./(1 + 4/3*pi*R.^3.*Nb_0);
Nb = Nb_0./(1 + 4/3*pi*(R.^3-R_0.^3).*Nb_0);

%==========================================================================
%ODE solver listener functions
%==========================================================================
function [value,isterminal,direction] = eventDetection(t,X)
timef = toc;
%a. Define the timeout in seconds
%TimeOut = 1000000;
%timeFlag = toc-TimeOut;

%if X(end,:) >1e-14
%    sizeFlag = 0;
%else
%    sizeFlag=1;
%    disp('bubble too small');
%end


%if ((sizeFlag==1)||(timeFlag<=0))
%    value = 1;
%else
%    value = 0;
%end

if X(end)<1e-10
    disp('Bubble too small')
end

if timef>1e4
    disp('Timeout')
end

value = (X(end)-1e-10)*(timef-1e4);
isterminal = 1;
direction = 0;
