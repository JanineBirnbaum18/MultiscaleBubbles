function [DarcyFun,PermFun,WaterViscModel,OutgasFun] = getFunctions_outgas(Geometry,PermModel,OutgasModel)

switch Geometry
    case 'Radial'
        DarcyFun = @(m0_fun,M1,P,P0,radius,z_p,z_u,K,mu,rho,Nb,R,T,dt,...
            min_density)Spherical_perm(m0_fun,M1,P,P0,radius,z_p,z_u,K,mu,rho,...
            Nb,R,T,dt,min_density);
        switch OutgasModel
            case 'Diffusive'
                OutgasFun = @(H2O1, H2O2, K, z_T, dt1, dt2, BC, timescheme)Spherical_outgas(H2O1, H2O2, K, z_T, dt1, dt2, BC, timescheme);
            case 'None'
                OutgasFun = @(H2O1, H2O2, K, z_T, dt1, dt2, BC, timescheme)(H2O1);
        end
    case 'Cylindrical'
        DarcyFun = @(m0_fun,M1,P,P0,radius,z_p,z_u,K,mu,rho,Nb,R,T,dt,...
            min_density)Cylindrical_perm(m0_fun,M1,P,P0,radius,z_p,z_u,K,mu,...
            rho,Nb,R,T,dt,min_density);
        switch OutgasModel
            case 'Diffusive'
                OutgasFun = @(H2O1, H2O2, K, z_T, dt1, dt2, BC, timescheme)Cylindrical_outgas(H2O1, H2O2, K, z_T, dt1, dt2, BC, timescheme);
            case 'None'
                OutgasFun = @(H2O1, H2O2, K, z_T, dt1, dt2, BC, timescheme)(H2O1);
        end
end

switch PermModel
    case 'Mueller2005Eff'
        PermFun = @(phi,Cc)Mueller2005Eff(phi);

    case 'Mueller2005Exp'
        PermFun = @(phi,Cc)Mueller2005Exp(phi);
    case 'None'
        PermFun = @(phi,Cc)(0.*phi);
end

WaterViscModel = @(rho,T)IAPSViscModel(rho,T);

% permeable outgassing
function [m_loss] = Cylindrical_perm(m0_fun,M1,P,P0,radius,z_p,z_u,K,mu,rho,Nb,R,T,...
    dt,min_density)
[h1,h2,A,B,C,D,E,F] = FDcoeff([z_p,z_u(end)]);

ddz = diag(-D(2:end-1),-1) + diag(-E(1:end-1)) + diag(C(1:end-2),1);
ddz(1,1:3) = [-A(1), B(1), -C(1)];
ddz(end,end-1:end) = [-D(end-1), -E(end-1)];


PP = (dt*8.314.*T./(Nb.*(z_p - [0, z_p(1:end-1)]).*4/3*pi.*R.^3*18.015/1000).*K.*rho./mu)'.*ddz;
Pr = zeros(size(z_p));
Pr(end) = (dt*8.314.*T(end)./(Nb(end).*(z_p(end) -  z_p(end-1)).*4/3*pi.*R(end).^3*18.015/1000).*K(end).*rho(end)./mu(end))*C(end-1)*P0(end);

M = eye(size(PP)) - PP;
b = (P' + Pr');

P_loss = (M\b)';

m_loss = 0*z_p;
for j = 1:length(z_p)
    m_loss(j) = m0_fun(R(j),P_loss(j),T(j));
end
m_loss(m_loss<rho.*4/3*pi.*R.^3) = rho(m_loss<rho.*4/3*pi.*R.^3).*4/3*pi.*R(m_loss<rho.*4/3*pi.*R.^3).^3;
%if any(M<rho_gas.*4/3*pi.*R.^3)
%    i = 0;
%    while norm(M0-M)>1e-18 && i<10
%        M0 = M;

%        P = pb_fun(M,T,R);
%        M_temp = M1 - (dt*1./(Nb.*([z_p(2:end) z_u(end)] - [0, z_p(1:end-1)])).*K.*rho./mu)'.*dPdz;
%        M_temp(M_temp<rho_gas.*4/3*pi.*R.^3) = rho_gas(M_temp<rho_gas.*4/3*pi.*R.^3).*4/3*pi.*R(M<rho_gas.*4/3*pi.*R.^3).^3;
%        i = i + 1;
%    end
%else
%    M(M<min_density) = min_density(M<min_density);
%end

function [M] = Spherical_perm(m0_fun,M1,P,P0,radius,z_p,z_u,K,mu,rho,Nb,R,T,...
    dt,min_density)
[h1,h2,A,B,C,D,E,F] = FDcoeff([z_p,z_u(end)]);

dPdz = diag(-D(2:end-1),-1) + diag(-E(1:end-1)) + diag(C(1:end-2),1);
dPdz(1,1:3) = [-A(1), B(1), -C(1)];
dPdz(end,end-1:end) = [-D(end-1), -E(end-1)];

dPdz = dPdz.*((1000./18.015).*8.314.*T)./((4.*pi./3).*R.^3)';
dPdz_r = C(end-1).*((1000./18.015).*8.314.*T)./((4.*pi./3).*R.^3);
dPdz_r(1:end) = 0;

MM = eye(size(dPdz)) - (dt*3*z_p.^2./(Nb.*(z_p.^3 - [0, z_p(1:end-1)].^3)).*K.*rho./mu)'.*dPdz;
Mr = 3*z_p.^2./(Nb.*(z_p.^3 - [0, z_p(1:end-1)].^3)).*K.*rho./mu.*dPdz_r;
M = ((MM)\(M1 + dt*Mr)')';
M0 = 0*M;

if any(M>min_density)
    i = 0;
    while norm(M0-M)>1e-18 && i<10
        M0 = M;
    
        P = pb_fun(M,T,R);
        dPdz = [-A(1)*P(1) + B(1)*P(2) - C(1)*P(3), ...
            -D(2:end-1).*P(1:end-2) - E(2:end-1).*P(2:end-1) + C(2:end-1).*P(3:end), ...
            -D(end).*P(end-1) - E(end).*P(end) + C(end).*P0(end)];
    
        M_temp = M1 + dt*3*z_p.^2./(Nb.*(z_p.^3 - [0, z_p(1:end-1)].^3)).*K.*rho./mu.*dPdz;
        M(M0>min_density) = M_temp(M0>min_density);
        M(M<min_density) = min_density(M<min_density);
        i = i + 1;
    end
else
    M(M<min_density) = min_density(M<min_density);
end


% diffusive outgassing
function [H2O] = Cylindrical_outgas(H2O1, H2O2, K, z_T, dt1, dt2, BC, timescheme)

Bt = (dt1 + dt2)/dt1/dt2;
Dt = dt1/dt2/(dt1+dt2);
Ft = (dt2+2*dt1)/dt1/(dt2+dt1);

% [h1,h2,A,B,C,D,E,F] = FDcoeff(z_T);
% dH2Odz = diag(-D(2:end),-1) + diag(-E) + diag(C(1:end-1),1);
% dH2Odz(1,1:3) = [-A(1), B(1), -C(1)];
% dH2Odz(end,end-2:end) = [D(end), -B(end), F(end)];
% 
% dKdz = [-A(1)*K(1) + B(1)*K(2) - C(1)*K(3), ...
%         -D(2:end-1).*K(1:end-2) - E(2:end-1).*K(2:end-1) + C(2:end-1).*K(3:end), ...
%         D(end).*K(end-2) - B(end).*K(end-1) + F(end).*K(end)];
% 
% d2H2Odz2 = diag(2*h2(2:end)./(h1(2:end).*h2(2:end).*(h1(2:end)+h2(2:end))),-1) +...
%     diag(-2*(h1+h2)./(h1.*h2.*(h1+h2))) + ...
%     diag(2*h1(1:end-1)./(h1(1:end-1).*h2(1:end-1).*(h1(1:end-1)+h2(1:end-1))),1);
% d2H2Odz2(1,1:4) = [2*(3*h1(1) + 2*h2(1) + h2(3))./h1(1)./(h1(1) + h2(1))./(h1(1)+h2(1)+h2(3)), ...
%                  -2*(2*h1(1) + 2*h2(1) + h2(3))./h1(1)./h2(1)./(h2(1) + h2(3)),...
%                  2*(2*h1(1) + h2(1) + h2(3))./(h1(1) + h2(1))./h2(1)./h2(3),...
%                  -2*(2*h1(1) + h2(1))./(h1(1) + h2(1) + h2(3))./(h2(1) + h2(3))./h2(3)];
% d2H2Odz2(end,end-3:end) = [-2*(h2(end-2) + 2*h2(end))./h1(end-2)./(h1(end-2)+h2(end-2))./(h1(end-2) + h2(end-2) + h2(end)),...
%                          2*(h1(end-2) + h2(end-2) + 2*h2(end))./h1(end-2)./h2(end-2)./(h2(end-2)+h2(end)),...
%                          -2*(h1(end-2) + 2*h2(end-2) + 2*h2(end))./(h1(end-2)+h2(end-2))./h2(end-2)./h2(end),...
%                          2*(h1(end-2) + 2*h2(end-2) + 3*h2(end))./(h1(end-2) + h2(end-2) + h2(end))./(h2(end-2) + h2(end))./h2(end)];


%XX = (dKdz'.*dH2Odz + K'.*d2H2Odz2);

Kstag = (K(2:end) + K(1:end-1))/2;
zstag = (z_T(2:end) + z_T(1:end-1))/2;
dH2Odt = -[0 1./diff(zstag) 0]'.*(diag([-Kstag(1:end-1)./(z_T(2:end-1) - z_T(1:end-2)), 0],-1) + ...
    diag([0, Kstag(1:end-1)./(z_T(2:end-1) - z_T(1:end-2)), 0]) + ...
    + diag([0, Kstag(2:end)./(z_T(3:end) - z_T(2:end-1)), 0]) + ...
    - diag([0, Kstag(2:end)./(z_T(3:end) - z_T(2:end-1))],1));
dH2Odt(1,:) = dH2Odt(2,:);
dH2Odt(end,:) = dH2Odt(end-1,:);

switch timescheme
    case 'BDF1'
    M = eye(size(dH2Odt)) - dt1.*dH2Odt;
    b = H2O1';

    case 'BDF2'
        M = eye(size(dH2Odt)) - 1/Ft.*dH2Odt;
        b = (-Dt/Ft.*H2O2' + Bt/Ft.*H2O1');

    case 'Steady'
        M = -dH2Odt;
        b = 0*H2O1';
end


% Dirichlet boundary
M(end,1:end) = 0;
M(end,end) = 1;
b(end) = BC;

% Symmetry boundary
M(1,:) = 0;
M(1,1:2) = [1, -1];
b(1) = 0;

H2O = (M\b)';

function [H2O] = Spherical_outgas(H2O1, H2O2, K, z_T, dt1, dt2, BC, timescheme)
Bt = (dt1 + dt2)/dt1/dt2;
Dt = dt1/dt2/(dt1+dt2);
Ft = (dt2+2*dt1)/dt1/(dt2+dt1);

Kstag = (K(2:end) + K(1:end-1))/2;
zstag = (z_T(2:end) + z_T(1:end-1))/2;
dH2Odt = -[0 1./(z_T(2:end-1).^2)./diff(zstag) 0]'.*(diag([-zstag(1:end-1).^2.*Kstag(1:end-1)./(z_T(2:end-1) - z_T(1:end-2)), 0],-1) + ...
    diag([0, zstag(1:end-1).^2.*Kstag(1:end-1)./(z_T(2:end-1) - z_T(1:end-2)), 0]) + ...
    + diag([0, zstag(2:end).^2.*Kstag(2:end)./(z_T(3:end) - z_T(2:end-1)), 0]) + ...
    - diag([0, zstag(2:end).^2.*Kstag(2:end)./(z_T(3:end) - z_T(2:end-1))],1));
dH2Odt(1,:) = dH2Odt(2,:);
dH2Odt(end,:) = dH2Odt(end-1,:);

%dH2Odt = [dH2Odt(1) dH2Odt 0];


% [h1,h2,A,B,C,D,E,F] = FDcoeff(z_T);
% dH2Odz = diag(-D(2:end),-1) + diag(-E) + diag(C(1:end-1),1);
% dH2Odz(1,1:3) = [-A(1), B(1), -C(1)];
% dH2Odz(end,end-2:end) = [D(end), -B(end), F(end)];
% 
% dKdz = [-A(1)*K(1) + B(1)*K(2) - C(1)*K(3), ...
%         -D(2:end-1).*K(1:end-2) - E(2:end-1).*K(2:end-1) + C(2:end-1).*K(3:end), ...
%         D(end).*K(end-2) - B(end).*K(end-1) + F(end).*K(end)];
% 
% d2H2Odz2 = diag(2*h2(2:end)./(h1(2:end).*h2(2:end).*(h1(2:end)+h2(2:end))),-1) +...
%     diag(-2*(h1+h2)./(h1.*h2.*(h1+h2))) + ...
%     diag(2*h1(1:end-1)./(h1(1:end-1).*h2(1:end-1).*(h1(1:end-1)+h2(1:end-1))),1);
% d2H2Odz2(1,1:4) = [2*(3*h1(1) + 2*h2(1) + h2(3))./h1(1)./(h1(1) + h2(1))./(h1(1)+h2(1)+h2(3)), ...
%                  -2*(2*h1(1) + 2*h2(1) + h2(3))./h1(1)./h2(1)./(h2(1) + h2(3)),...
%                  2*(2*h1(1) + h2(1) + h2(3))./(h1(1) + h2(1))./h2(1)./h2(3),...
%                  -2*(2*h1(1) + h2(1))./(h1(1) + h2(1) + h2(3))./(h2(1) + h2(3))./h2(3)];
% d2H2Odz2(end,end-3:end) = [-2*(h2(end-2) + 2*h2(end))./h1(end-2)./(h1(end-2)+h2(end-2))./(h1(end-2) + h2(end-2) + h2(end)),...
%                          2*(h1(end-2) + h2(end-2) + 2*h2(end))./h1(end-2)./h2(end-2)./(h2(end-2)+h2(end)),...
%                          -2*(h1(end-2) + 2*h2(end-2) + 2*h2(end))./(h1(end-2)+h2(end-2))./h2(end-2)./h2(end),...
%                          2*(h1(end-2) + 2*h2(end-2) + 3*h2(end))./(h1(end-2) + h2(end-2) + h2(end))./(h2(end-2) + h2(end))./h2(end)];


%%XX = ((2*K./z_T + dKdz)'.*dH2Odz + K'.*d2H2Odz2);
%XX = dH2Odz*(z_T.^2.*K'.*dH2Odz)./z_T.^2;
switch timescheme
    case 'BDF1'
    M = eye(size(dH2Odt)) - dt1.*dH2Odt;
    b = H2O1';

    case 'BDF2'
        M = eye(size(dH2Odt)) - 1/Ft.*dH2Odt;
        b = (-Dt/Ft.*H2O2' + Bt/Ft.*H2O1');

    case 'Steady'
        M = -dH2Odt;
        b = 0*H2O1';   
end

% Dirichlet boundary
M(end,1:end) = 0;
M(end,end) = 1;
b(end) = BC;

% Symmetry boundary
M(1,:) = 0;
M(1,1:2) = [-1,1];
b(1) = 0;

H2O = (M\b)';
    
% Permeability
function [k] = Mueller2005Eff(phi)
k = 1e-17*((100*phi).^3.4);

function [k] = Mueller2005Exp(phi)
k = [1e-15*(100*phi-30).^2];
k(phi<=0.3001) = 1e-25;

% Water vapor viscosity
function [mu] = IAPSViscModel(rho,T)
Hi = [1.67752, 2.20462, 0.6366564, -0.241605];
T0 = 647.096;
rho0 = 322.0; 

mu0 = 0;
for i = 0:3
    mu0 = mu0 + (Hi(i+1)./(T/T0).^i);
end
mu0 = 100*sqrt(T/T0)./mu0;

Hij = zeros(6,7);
Hij(1,1) = 5.20094e-1; 
Hij(2,1) = 8.50895e-2;
Hij(3,1) = -1.08374;
Hij(4,1) = -2.89555e-1;
Hij(1,2) = 2.22531e-1;
Hij(2,2) = 9.99115e-1;
Hij(3,2) = 1.88797;
Hij(4,2) = 1.26613;
Hij(6,2) = 1.20573e-1;
Hij(1,3) = -2.81378e-1;
Hij(2,3) = -9.06851e-1;
Hij(3,3) = -7.724792e-1;
Hij(4,3) = -4.89837e-1;
Hij(5,3) = -2.57040e-1;
Hij(1,4) = 1.61913e-1;
Hij(2,4) = 2.57399e-1;
Hij(1,5) = -3.25372e-2;
Hij(4,5) = 6.98452e-2;
Hij(5,6) = 8.72102e-3;
Hij(4,7) = -4.35673e-3;
Hij(6,7) = -5.93264e-4;

mu1 = 0;
for i = 0:5
    for j = 0:6
        mu1 = mu1 + (1./(T/T0)-1).^(i).*Hij(i+1,j+1).*(rho/rho0-1).^j;
    end
end
mu = mu0.*exp(rho./rho0.*mu1).*1e-6;

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
