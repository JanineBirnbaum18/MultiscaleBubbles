function [ThermFun,cpmeltFun,cpFun,kmeltFun,kFun,rhoFun] = getFunctions_thermal(Geometry,BC_type,CpModel,kmodel,rhoModel,alpha,T0)

switch Geometry
    case 'Radial'
        ThermFun = @(T1,T2,rho,cp,k,z_t,dt1,dt2,BC,flux,timescheme)Spherical_temp(T1,T2,rho,cp,k,z_t,dt1,dt2,BC_type,BC,flux,timescheme);

    case 'Cylindrical'
        ThermFun = @(T1,T2,rho,cp,k,z_t,dt1,dt2,BC,flux,timescheme)Cylindrical_temp(T1,T2,rho,cp,k,z_t,dt1,dt2,BC_type,BC,flux,timescheme);
   
    case 'Cylindrical-2D'
        ThermFun = @(T1,T2,rho,cp,k,z_t,r_t,dt1,dt2,BC,flux,timescheme)Cylindrical2D_temp(T1,T2,rho,cp,k,z_t,r_t,dt1,dt2,BC_type,BC,flux,timescheme);
end

switch CpModel
    case 'BagdassarovDingwell1994'
        cpFun = @(phi,cpmelt,rhomelt,T,P)BagdassarovDingwell1994Cp(phi,cpmelt,rhomelt,T,P);
end

switch rhoModel
    case 'BagdassarovDingwell1994'
        rhoFun = @(rho0,T)BagdassarovDingwell1994rho(rho0,T);
    case 'ConstantExpansivity'
        rhoFun = @(rho0,T)ConstantExpansivity(rho0,T,alpha,T0);
end

kmeltFun = @(Composition)BagdassarovDingwell1994kmelt(Composition);
kFun = @(phi,kmelt)BagdassarovDingwell1994kfoam(phi,kmelt,kmodel);
rhoFun = @(rho0,T)BagdassarovDingwell1994rho(rho0,T);
cpmeltFun = @(Composition,T,wtH2O)Stebbins(Composition,T,wtH2O);

function [T] = Cylindrical_temp(T1,T2,rho,cp,k,z_T,dt1,dt2,BC_type,BC,flux,timescheme)

Bt = (dt1 + dt2)/dt1/dt2;
Dt = dt1/dt2/(dt1+dt2);
Ft = (dt2+2*dt1)/dt1/(dt2+dt1);

n = 1;
tol = 1e-5;
max_iter = 10;
T = T1;
Ti = T1;

[h1,h2,A,B,C,D,E,F] = FDcoeff(z_T);
dTdz = diag(-D(2:end),-1) + diag(-E) + diag(C(1:end-1),1);
dTdz(1,1:3) = [-A(1), B(1), -C(1)];
dTdz(end,end-2:end) = [D(end), -B(end), F(end)];

dkdz = [-A(1)*k(1) + B(1)*k(2) - C(1)*k(3), ...
        -D(2:end-1).*k(1:end-2) - E(2:end-1).*k(2:end-1) + C(2:end-1).*k(3:end), ...
        D(end).*k(end-2) - B(end).*k(end-1) + F(end).*k(end)];

d2Tdz2 = diag(2*h2(2:end)./(h1(2:end).*h2(2:end).*(h1(2:end)+h2(2:end))),-1) +...
    diag(-2*(h1+h2)./(h1.*h2.*(h1+h2))) + ...
    diag(2*h1(1:end-1)./(h1(1:end-1).*h2(1:end-1).*(h1(1:end-1)+h2(1:end-1))),1);
d2Tdz2(1,1:4) = [2*(3*h1(1) + 2*h2(1) + h2(3))./h1(1)./(h1(1) + h2(1))./(h1(1)+h2(1)+h2(3)), ...
             -2*(2*h1(1) + 2*h2(1) + h2(3))./h1(1)./h2(1)./(h2(1) + h2(3)),...
             2*(2*h1(1) + h2(1) + h2(3))./(h1(1) + h2(1))./h2(1)./h2(3),...
             -2*(2*h1(1) + h2(1))./(h1(1) + h2(1) + h2(3))./(h2(1) + h2(3))./h2(3)];
d2Tdz2(end,end-3:end) = [-2*(h2(end-2) + 2*h2(end))./h1(end-2)./(h1(end-2)+h2(end-2))./(h1(end-2) + h2(end-2) + h2(end)),...
                     2*(h1(end-2) + h2(end-2) + 2*h2(end))./h1(end-2)./h2(end-2)./(h2(end-2)+h2(end)),...
                     -2*(h1(end-2) + 2*h2(end-2) + 2*h2(end))./(h1(end-2)+h2(end-2))./h2(end-2)./h2(end),...
                     2*(h1(end-2) + 2*h2(end-2) + 3*h2(end))./(h1(end-2) + h2(end-2) + h2(end))./(h2(end-2) + h2(end))./h2(end)];

TT = (1./rho./cp)'.*(dkdz'.*dTdz + k'.*d2Tdz2);
TR = (1./rho./cp)'.*flux;

switch timescheme
    case 'BDF1'
    M = eye(size(TT)) - dt1.*TT;
    b = T1' + dt1.*TR;

    case 'BDF2'
        M = eye(size(TT)) - 1/Ft.*TT;
        b = (-Dt/Ft.*T2' + Bt/Ft.*T1' + 1/Ft.*TR);

    case 'Steady'
        M = -TT;
        b = 0*T1';
end

% No flux at bottom
M(1,:) = 0;
M(1,1:3) = [-A(1), B(1), -C(1)];
b(1) = 0;

while (n<max_iter && norm(T-Ti)>tol) || n==1
    
    % Applied temperature at top
    switch BC_type
        case 'Dirichlet'
    
        M(end,:) = 0;
        M(end,end) = 1;
        b(end) = BC;
    
        case 'Forced'
        M(end,:) = 0;
        M(end,end-2:end) = k(end).*[D(end), -B(end), F(end)] + [0,0,BC(1) + 0*5.67e-8.*Ti(end).^3];
        b(end) = BC(1).*BC(2) + 0*5.67e-8.*BC(2).^4;
    
    end
    
    Ti = T;
    T = (M\b)';
    
    n = n+1;
    err = norm(T-Ti);

end

function [T] = Spherical_temp(T1,T2,rho,cp,k,z_T,dt1,dt2,BC_type,BC,flux,timescheme)

%[~,~,~,Bt,~,Dt,~,Ft] = FDcoeff([0,dt2,dt1+dt2]);
Bt = (dt1 + dt2)/dt1/dt2;
Dt = dt1/dt2/(dt1+dt2);
Ft = (dt2+2*dt1)/dt1/(dt2+dt1);

n = 1;
tol = 1e-5;
max_iter = 10;
T = T1;
Ti = T1;

[h1,h2,A,B,C,D,E,F] = FDcoeff(z_T);
dTdz = diag(-D(2:end),-1) + diag(-E) + diag(C(1:end-1),1);
dTdz(1,1:3) = [-A(1), B(1), -C(1)];
dTdz(end,end-2:end) = [D(end), -B(end), F(end)];

dkdz = [-A(1)*k(1) + B(1)*k(2) - C(1)*k(3), ...
        -D(2:end-1).*k(1:end-2) - E(2:end-1).*k(2:end-1) + C(2:end-1).*k(3:end), ...
        D(end).*k(end-2) - B(end).*k(end-1) + F(end).*k(end)];

d2Tdz2 = diag(2*h2(2:end)./(h1(2:end).*h2(2:end).*(h1(2:end)+h2(2:end))),-1) +...
    diag(-2*(h1+h2)./(h1.*h2.*(h1+h2))) + ...
    diag(2*h1(1:end-1)./(h1(1:end-1).*h2(1:end-1).*(h1(1:end-1)+h2(1:end-1))),1);
d2Tdz2(1,1:4) = [2*(3*h1(1) + 2*h2(1) + h2(3))./h1(1)./(h1(1) + h2(1))./(h1(1)+h2(1)+h2(3)), ...
             -2*(2*h1(1) + 2*h2(1) + h2(3))./h1(1)./h2(1)./(h2(1) + h2(3)),...
             2*(2*h1(1) + h2(1) + h2(3))./(h1(1) + h2(1))./h2(1)./h2(3),...
             -2*(2*h1(1) + h2(1))./(h1(1) + h2(1) + h2(3))./(h2(1) + h2(3))./h2(3)];
d2Tdz2(end,end-3:end) = [-2*(h2(end-2) + 2*h2(end))./h1(end-2)./(h1(end-2)+h2(end-2))./(h1(end-2) + h2(end-2) + h2(end)),...
                     2*(h1(end-2) + h2(end-2) + 2*h2(end))./h1(end-2)./h2(end-2)./(h2(end-2)+h2(end)),...
                     -2*(h1(end-2) + 2*h2(end-2) + 2*h2(end))./(h1(end-2)+h2(end-2))./h2(end-2)./h2(end),...
                     2*(h1(end-2) + 2*h2(end-2) + 3*h2(end))./(h1(end-2) + h2(end-2) + h2(end))./(h2(end-2) + h2(end))./h2(end)];


TT = (1./rho./cp)'.*((2*k./z_T + dkdz)'.*dTdz + k'.*d2Tdz2);

switch timescheme
    case 'BDF1'
    M = eye(size(TT)) - dt1.*TT;
    b = T1';

    case 'BDF2'
        M = eye(size(TT)) - 1/Ft.*TT;
        b = (-Dt/Ft.*T2' + Bt/Ft.*T1');

    case 'Steady'
        M = -TT;
        b = 0*T1';
end

% Symmetry in temperature
M(1,:) = 0;
M(1,1:3) = [-A(1), B(1), -C(1)];
b(1) = 0;

while (n<max_iter && norm(T-Ti)>tol) || n==1
    % Applied temperature at outer edge
    switch BC_type
        case 'Dirichlet'
    
        M(end,1:end) = 0;
        M(end,end) = 1;
        b(end) = BC;
    
        case 'Forced'
        M(end,:) = 0;
        M(end,end-2:end) = k(end).*[D(end), -B(end), F(end)] + [0,0,BC(1) + 0*5.67e-8.*Ti(end-1).^3];
        b(end) = BC(1).*BC(2) + 0*5.67e-8.*BC(2).^4;
    end
    
    Ti = T;
    T = (M\b)';
    
    n = n+1;
    err = norm(T-Ti);
end

function [T] = Cylindrical2D_temp(T1,T2,rho,cp,k,z_T,r_T,dt1,dt2,BC_type,BC,flux,timescheme)

Bt = (dt1 + dt2)/dt1/dt2;
Dt = dt1/dt2/(dt1+dt2);
Ft = (dt2+2*dt1)/dt1/(dt2+dt1);

n = 1;
tol = 1e-5;
max_iter = 10;
T = T1;
Ti = T1;

while (n<max_iter && norm(T-Ti)>tol) || n==1
 
    [h1,h2,A,B,C,D,E,F] = FDcoeff(r_T);
    dTdr = diag(-D(2:end),-1) + diag(-E) + diag(C(1:end-1),1);
    dTdr(1,1:3) = [-A(1), B(1), -C(1)];
    dTdr(end,end-2:end) = [D(end), -B(end), F(end)];
    
    dkdr = [-A(1)*k(1) + B(1)*k(2) - C(1)*k(3), ...
            -D(2:end-1).*k(1:end-2) - E(2:end-1).*k(2:end-1) + C(2:end-1).*k(3:end), ...
            D(end).*k(end-2) - B(end).*k(end-1) + F(end).*k(end)];
    
    d2Tdr2 = diag(2*h2(2:end)./(h1(2:end).*h2(2:end).*(h1(2:end)+h2(2:end))),-1) +...
        diag(-2*(h1+h2)./(h1.*h2.*(h1+h2))) + ...
        diag(2*h1(1:end-1)./(h1(1:end-1).*h2(1:end-1).*(h1(1:end-1)+h2(1:end-1))),1);
    d2Tdr2(1,1:4) = [2*(3*h1(1) + 2*h2(1) + h2(3))./h1(1)./(h1(1) + h2(1))./(h1(1)+h2(1)+h2(3)), ...
                 -2*(2*h1(1) + 2*h2(1) + h2(3))./h1(1)./h2(1)./(h2(1) + h2(3)),...
                 2*(2*h1(1) + h2(1) + h2(3))./(h1(1) + h2(1))./h2(1)./h2(3),...
                 -2*(2*h1(1) + h2(1))./(h1(1) + h2(1) + h2(3))./(h2(1) + h2(3))./h2(3)];
    d2Tdr2(end,end-3:end) = [-2*(h2(end-2) + 2*h2(end))./h1(end-2)./(h1(end-2)+h2(end-2))./(h1(end-2) + h2(end-2) + h2(end)),...
                         2*(h1(end-2) + h2(end-2) + 2*h2(end))./h1(end-2)./h2(end-2)./(h2(end-2)+h2(end)),...
                         -2*(h1(end-2) + 2*h2(end-2) + 2*h2(end))./(h1(end-2)+h2(end-2))./h2(end-2)./h2(end),...
                         2*(h1(end-2) + 2*h2(end-2) + 3*h2(end))./(h1(end-2) + h2(end-2) + h2(end))./(h2(end-2) + h2(end))./h2(end)];

    [h1,h2,A,B,C,D,E,F] = FDcoeff(z_T);
    dTdz = diag(-D(2:end),-1) + diag(-E) + diag(C(1:end-1),1);
    dTdz(1,1:3) = [-A(1), B(1), -C(1)];
    dTdz(end,end-2:end) = [D(end), -B(end), F(end)];
    
    dkdz = [-A(1)*k(1) + B(1)*k(2) - C(1)*k(3), ...
            -D(2:end-1).*k(1:end-2) - E(2:end-1).*k(2:end-1) + C(2:end-1).*k(3:end), ...
            D(end).*k(end-2) - B(end).*k(end-1) + F(end).*k(end)];
    
    d2Tdz2 = diag(2*h2(2:end)./(h1(2:end).*h2(2:end).*(h1(2:end)+h2(2:end))),-1) +...
        diag(-2*(h1+h2)./(h1.*h2.*(h1+h2))) + ...
        diag(2*h1(1:end-1)./(h1(1:end-1).*h2(1:end-1).*(h1(1:end-1)+h2(1:end-1))),1);
    d2Tdz2(1,1:4) = [2*(3*h1(1) + 2*h2(1) + h2(3))./h1(1)./(h1(1) + h2(1))./(h1(1)+h2(1)+h2(3)), ...
                 -2*(2*h1(1) + 2*h2(1) + h2(3))./h1(1)./h2(1)./(h2(1) + h2(3)),...
                 2*(2*h1(1) + h2(1) + h2(3))./(h1(1) + h2(1))./h2(1)./h2(3),...
                 -2*(2*h1(1) + h2(1))./(h1(1) + h2(1) + h2(3))./(h2(1) + h2(3))./h2(3)];
    d2Tdz2(end,end-3:end) = [-2*(h2(end-2) + 2*h2(end))./h1(end-2)./(h1(end-2)+h2(end-2))./(h1(end-2) + h2(end-2) + h2(end)),...
                         2*(h1(end-2) + h2(end-2) + 2*h2(end))./h1(end-2)./h2(end-2)./(h2(end-2)+h2(end)),...
                         -2*(h1(end-2) + 2*h2(end-2) + 2*h2(end))./(h1(end-2)+h2(end-2))./h2(end-2)./h2(end),...
                         2*(h1(end-2) + 2*h2(end-2) + 3*h2(end))./(h1(end-2) + h2(end-2) + h2(end))./(h2(end-2) + h2(end))./h2(end)];


    
    TT = (1./rho./cp)'.*(dkdz'.*dTdz + k'.*d2Tdz2);
    TR = (1./rho./cp)'.*flux;
    
    switch timescheme
        case 'BDF1'
        M = eye(size(TT)) - dt1.*TT;
        b = T1' + dt1.*TR;
    
        case 'BDF2'
            M = eye(size(TT)) - 1/Ft.*TT;
            b = (-Dt/Ft.*T2' + Bt/Ft.*T1' + 1/Ft.*TR);
    
        case 'Steady'
            M = -TT;
            b = 0*T1';
    end
    
    % Applied temperature at top
    switch BC_type
        case 'Dirichlet'
    
        M(end,:) = 0;
        M(end,end) = 1;
        b(end) = BC;
    
        case 'Forced'
        M(end,:) = 0;
        M(end,end-2:end) = k(end).*[D(end), -B(end), F(end)] + [0,0,BC(1) + 0*5.67e-8.*Ti(end).^3];
        b(end) = BC(1).*BC(2) + 0*5.67e-8.*BC(2).^4;
    
    end
    
    % No flux at bottom
    M(1,:) = 0;
    M(1,1:3) = [-A(1), B(1), -C(1)];
    b(1) = 0;
    
    Ti = T;
    T = (M\b)';
    
    n = n+1;
    err = norm(T-Ti);

end

% finite differences coefficients
function [h1,h2,A,B,C,D,E,F] = FDcoeff(z)
h1 = [z(2)-z(1), z(2:end-1)-z(1:end-2), z(end-1)-z(end-2)];
h2 = [z(3)-z(2), z(3:end)-z(2:end-1), z(end)-z(end-1)];
A = (2*h1 + h2)./h1./(h1+h2);
B = (h1+h2)./h1./h2;
C = h1./(h1+h2)./h2;
D = h2./h1./(h1+h2);
E = (h1-h2)./h1./h2; 
F = (h1 + 2*h2)./h2./(h1+h2);


% thermal material properties
function [cpmelt] = Stebbins(Composition,T,wtH2O)

Composition = Composition'./(sum(Composition) + wtH2O/100);
Composition(11,:) = wtH2O'/100;

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
mH = 1.00784;
mF = 19.00;
mO = 15.999;

%Molar mass (g/mol) of oxides
OxideMolarMass = zeros(12,1);
OxideMolarMass(1) =  (mSi+2*mO);
OxideMolarMass(2) = (mTi+2*mO);
OxideMolarMass(3) = (2*mAl+3*mO);
OxideMolarMass(4) = (1*mFe+1*mO);
OxideMolarMass(5) = (1*mMn+1*mO);
OxideMolarMass(6) = (1*mMg+1*mO);
OxideMolarMass(7) = (1*mCa+1*mO);
OxideMolarMass(8) = (2*mNa+1*mO);
OxideMolarMass(9) = (2*mK+1*mO);
OxideMolarMass(10) = (2*mP+5*mO);
OxideMolarMass(11) = (2*mH+mO);
OxideMolarMass(12) = (2*mF+mO);

Xhyd = Composition./OxideMolarMass;
XH2O = Xhyd(11,:)./sum(Xhyd,1);
Xanhyd = Xhyd;
Xanhyd(5,:) = 0;
Xanhyd(10:12,:) = 0;
Xanhyd = Xanhyd./sum(Xanhyd,1);

a = [66.354, 33.851, 91.404, 40.949, 0, 32.244, 46.677, 69.067, 107.194, 0, 0, 0]';
b = 1e-2*[0.7797, 6.4572, 4.4940, 2.9349, 0, 2.7288, 0.3565, 1.8603, -3.2194, 0, 0, 0]';
c = 1e5*[-28.003, 4.470, -21.465, -7.6986, 0, 1.7549, -1.9322, 2.9101, -28.929, 0, 0, 0]';

Cpanhyd = sum(a.*Xanhyd,1) + sum(b.*Xanhyd,1).*T + sum(c.*Xanhyd,1).*(T.^(-2));

M = sum(Xhyd./sum(Xhyd,1).*OxideMolarMass,1);
cpmelt = (XH2O.*(237) + (1-XH2O).*Cpanhyd).*M;

function [cpfoam] = BagdassarovDingwell1994Cp(phi,cpmelt,rhomelt,T,P)
cpgas = (-6e-08.*T.^2 + 0.0008.*T + 1.558)*1000;

coef = zeros(10,6);
coef(1,3)=0.24657688e6;
coef(1,4)=0.51359951e2;
coef(2,3)=0.58638965e0;
coef(2,4)=-0.28646939e-2;
coef(2,5)=0.31375577e-4;
coef(3,3)=-0.62783840e1;
coef(3,4)=0.14791599e-1;
coef(3,5)=0.35779579e-3;
coef(3,6)=0.15432925e-7;
coef(4,4)=-0.42719875e0;

coef(4,5)=-0.16325155e-4;
coef(5,3)=0.56654978e4;
coef(5,4)=-0.16580167e2;
coef(5,5)=0.76560762e-1;
coef(6,4)=0.10917883e0;
coef(7,1)=0.38878656e13;
coef(7,2)=-0.13494878e9;
coef(7,3)=0.30916564e6;
coef(7,4)=0.75591105e1;
coef(8,3)=-0.65537898e5;
coef(8,4)=0.18810675e3;
coef(9,1)=-0.14182435e14;
coef(9,2)=0.18165390e9;
coef(9,3)=-0.19769068e6;
coef(9,4)=-0.23530318e2;
coef(10,3)=0.92093375e5;
coef(10,4)=0.12246777e3;

P = P.*1e-5; %1 pascal a 1e-5 bars

coef2=zeros(10,length(T));
for i=1:10
    coef2(i,:)=coef(i,1)*T.^-4 + coef(i,2)*T.^-2 + coef(i,3)*T.^-1 +...
        coef(i,4) + coef(i,5)*T + coef(i,6)*T.^2;
end
% PRT = P/RT where P [bars], R [cm^3*bar/K/mol] and T [K]
PRT = P./(83.14472*T);
% solve implicit equation for rho and conver
% t to kg/m^3
fun = @(rho)PS_myfun(rho,coef2,PRT);
options = optimoptions('fsolve','Display','none','Algorithm','levenberg-marquardt');
rhogas = fsolve(fun,0.001*ones(size(T)),options)*18.01528*1000;

cpfoam = (phi.*cpgas.*rhogas + (1-phi).*cpmelt.*rhomelt)./rhomelt;

function y = PS_myfun(rho,coef2,PRT)
y = (rho+coef2(1,:).*rho.^2-rho.^2.*((coef2(3,:)+2*coef2(4,:).*rho+...
    3*coef2(5,:).*rho.^2+4*coef2(6,:).*rho.^3)./((coef2(2,:)+coef2(3,:).*rho+...
    coef2(4,:).*rho.^2+coef2(5,:).*rho.^3+coef2(6,:).*rho.^4).^2))+...
    coef2(7,:).*rho.^2.*exp(-coef2(8,:).*rho)+coef2(9,:).*rho.^2.*exp(-coef2(10,:).*rho)) - PRT;

function [kmelt] = BagdassarovDingwell1994kmelt(Composition)
%[SiO2 TiO2 Al2O3 FeO(T) MnO MgO CaO Na2O K2O P2O5 H2O F2O-1]
Pi = [-1.062, 0, 0.449, -1.687, 0, -3.318, -1.905, -1.952, -3.939, 0, 0, 0];
Wi = Composition./(sum(Composition([1,3, 4, 6, 7, 8, 9])));
kmelt = 2.371 + sum(Pi.*Wi);

function [k] = BagdassarovDingwell1994kfoam(phi,kmelt,model)

switch model
    case 'RayleighMaxwell'
    k = kmelt.*(1-phi)./(1+phi);
    case 'Bruggeman'
    k = kmelt.*(1-phi).^(3/2);
end

function [rho] = BagdassarovDingwell1994rho(rho0,T)
rho=rho0*(1 + 4.635e-6.*(T-273) + 0.654e-9.*(T-273).^2).^(-3); 

function[rho] = ConstantExpansivity(rho0,T,alpha, T0)
rho = rho0*(1 - (T-T0)*alpha);
