function [DynFun] = getFunctions_dynamic_dimensionless(Geometry,BC)

%warning('off','MATLAB:illConditionedMatrix')
%warning('off','MATLAB:nearlySingularMatrix')

switch Geometry
    case 'Radial'
        DynFun = @(P1,u1,P2,u2,P0,rho,drhodt,phi,a,SurfTens,eta,beta,dbetadt,...
            R,g,z_p,z_u,dt1,dt2,timescheme)Spherical_pde(P1,u1,P2,u2,P0,rho,drhodt,...
            phi,a,SurfTens,eta,beta,dbetadt,R,g,z_p,z_u,dt1,dt2,BC,timescheme);

    case 'Cylindrical'
        DynFun = @(P1,u1,P2,u2,P0,rho,drhodt,phi,a,SurfTens,eta,beta,dbetadt,...
            R,g,z_p,z_u,dt1,dt2,timescheme)Cylindrical_pde(P1,u1,P2,u2,P0,rho,drhodt,...
            phi,a,SurfTens,eta,beta,dbetadt,R,g,z_p,z_u,dt1,dt2,BC,timescheme);

    case '2D Cylindrical'
        DynFun = @(P1,u1,P2,u2,P0,rho,drhodt,phi,a,SurfTens,eta,beta,dbetadt,...
            R,g,z_p,z_u,dt1,dt2,timescheme)Cylindrical2D_pde(P1,u1,P2,u2,P0,rho,drhodt,...
            phi,a,SurfTens,eta,beta,dbetadt,R,g,z_p,z_u,dt1,dt2,BC,timescheme);

end

%%%%%%%%%%%%% Cylindrical %%%%%%%%%%%%%%%%%
function [P,u] = Cylindrical_pde(P1,u1,P2,u2,P0,rho,drhodt,...
            phi,a,SurfTens,eta,beta,dbetadt,R,g,z_p,z_u,dt1,dt2,BC,timescheme)

Mu = max(eta);
L = max(z_u);
Rho = max(rho);
U = max(1e-12,max(u1));
Re = U*L*Rho/Mu;

P1 = P1*(L/Mu/U);
u1 = u1/U;
P2 = P2*(L/Mu/U);
u2 = u2/U;
P0 = P0*(L/Mu/U);
rho = rho/Rho;
drhodt = drhodt/Rho*(L/U);
eta = eta/Mu;
beta = beta/(L/Mu/U);
z_p = z_p/L;
z_u = z_u/L;
dt1 = dt1/(L/U);
dt2 = dt2/(L/U);

%[~,~,~,Bt,~,Dt,~,Ft] = FDcoeff([0,dt2,dt1+dt2]);
Bt = (dt1 + dt2)/dt1/dt2;
Dt = dt1/dt2/(dt1+dt2);
Ft = (dt2+2*dt1)/dt1/(dt2+dt1);

P = P1; 
u = u1;
Vi = [P u]';
V = 0;

max_iter = 1000;
tol = 1e-7;
n = 1;

[h1,h2,A,B,C,D,E,F] = FDcoeff(z_p);
dPdz = diag(-D(2:end),-1) + diag(-E) + diag(C(1:end-1),1);
dPdz(1,1:3) = [-A(1), B(1), -C(1)];
dPdz(end,end-2:end) = [D(end), -B(end), F(end)];

dP0dz = [-A(1)*P0(1) + B(1)*P0(2) - C(1)*P0(3), ...
        -D(2:end-1).*P0(1:end-2) - E(2:end-1).*P0(2:end-1) + C(2:end-1).*P0(3:end), ...
        D(end).*P0(end-2) - B(end).*P0(end-1) + F(end).*P0(end)];
drhodz = [-A(1)*rho(2) + B(1)*rho(4) - C(1)*rho(6), ...
        -D(2:end-1).*rho(2:2:end-4) - E(2:end-1).*rho(4:2:end-2) + C(2:end-1).*rho(6:2:end), ...
        D(end).*rho(end-5) - B(end).*rho(end-3) + F(end).*rho(end-1)];
dbetadz = [-A(1)*beta(1) + B(1)*beta(2) - C(1)*beta(3), ...
        -D(2:end-1).*beta(1:end-2) - E(2:end-1).*beta(2:end-1) + C(2:end-1).*beta(3:end), ...
        D(end).*beta(end-2) - B(end).*beta(end-1) + F(end).*beta(end)];


h2 = z_u(2:end)-z_u(1:end-1);
dudz_P = diag(-1./[h2 h2(end)]) + diag(1./h2,1);
dudz_P = dudz_P(1:end-1,1:end);

h1 = z_u(2:end-1)-z_u(1:end-2);
h2 = z_u(3:end)-z_u(2:end-1);
w1 = (z_p - z_u(1:end-1))./[h1, h2(end)];
w2 = (z_u(2:end)-z_p)./[h1, h2(end)];
u_interp = diag([w1,w1(end)]) + diag(w2,1);
u_interp = u_interp(1:end-1,:);

[h1,h2,A,B,C,D,E,F] = FDcoeff(z_u);

dudz = diag(-D(2:end),-1) + diag(-E) + diag(C(1:end-1),1);
dudz(1,1:3) = [-A(1), B(1), -C(1)];
dudz(end,end-2:end) = [D(end), -B(end), F(end)];

h1 = z_u(2:end-1)-z_u(1:end-2);
h1 = [h1(1),h1,h1(end)];
h2 = z_u(3:end)-z_u(2:end-1);
h2 = [h2(1),h2,h2(end)];

d2udz2 = diag(2*h2(2:end)./(h1(2:end).*h2(2:end).*(h1(2:end)+h2(2:end))),-1) +...
    diag(-2*(h1+h2)./(h1.*h2.*(h1+h2))) + ...
    diag(2*h1(1:end-1)./(h1(1:end-1).*h2(1:end-1).*(h1(1:end-1)+h2(1:end-1))),1);
d2udz2(1,1:4) = [2*(3*h1(1) + 2*h2(1) + h2(3))./h1(1)./(h1(1) + h2(1))./(h1(1)+h2(1)+h2(3)), ...
                 -2*(2*h1(1) + 2*h2(1) + h2(3))./h1(1)./h2(1)./(h2(1) + h2(3)),...
                 2*(2*h1(1) + h2(1) + h2(3))./(h1(1) + h2(1))./h2(1)./h2(3),...
                 -2*(2*h1(1) + h2(1))./(h1(1) + h2(1) + h2(3))./(h2(1) + h2(3))./h2(3)];
d2udz2(end,end-3:end) = [-2*(h2(end-2) + 2*h2(end))./h1(end-2)./(h1(end-2)+h2(end-2))./(h1(end-2) + h2(end-2) + h2(end)),...
                         2*(h1(end-2) + h2(end-2) + 2*h2(end))./h1(end-2)./h2(end-2)./(h2(end-2)+h2(end)),...
                         -2*(h1(end-2) + 2*h2(end-2) + 2*h2(end))./(h1(end-2)+h2(end-2))./h2(end-2)./h2(end),...
                         2*(h1(end-2) + 2*h2(end-2) + 3*h2(end))./(h1(end-2) + h2(end-2) + h2(end))./(h2(end-2) + h2(end))./h2(end)];

h2 = [z_p(1) z_p(2:end) - z_p(1:end-1) z_u(end)-z_p(end)];
dPdz_u = diag(-1./h2(2:end),-1) + diag(1./h2);
dPdz_u(1,1:2) = dPdz_u(2,1:2);
dPdz_u = dPdz_u(:,1:end-1);

z_t = reshape([z_u; z_p, 0],1,[]);
z_t = z_t(1:end-1);

a_interp = griddedInterpolant(z_p,a,'linear','nearest');
a_interp = a_interp(z_t);
phi_interp = griddedInterpolant(z_p,phi,'linear','nearest');
phi_interp = phi_interp(z_t);

eta_interp = eta;
rho_interp = rho(1:2:end);

while (n<max_iter && norm(V-Vi)>tol) || n==1

PP = diag(-1./rho(2:2:end).*(drhodt + 1/2.*drhodz.*(u_interp*u')')) + ...
     diag(-1./beta.*(dbetadt + 1/2.*dbetadz.*(u_interp*u')')) + ...
     (-1/2*(u_interp*u').*(dPdz)) + ...
     diag(-1/2*(dudz_P*u')');
PR = (-1./rho(2:2:end).*(drhodt).*(1./beta-P0))' +...
    ((-1./beta).*dbetadt.*(-P0))';
PU = (-1./rho(2:2:end).*drhodz.*(1./beta + 1/2*P-P0))'.*u_interp + ...
    (-1./beta.*dbetadz.*(1/2*P-P0))'.*u_interp + ...
    (-1/2*(dPdz*P') + dP0dz').*u_interp + ...
    (-(1./beta + 1/2*P-P0))'.*dudz_P;

dudr = -3*u'./R;
dudr1 = -3*u1'./R;
dudr2 = -3*u2'./R;

switch timescheme
    case 'BDF1'
        d2udrdt = 0*dudr;
    case 'BDF2'
        d2udrdt = Ft*dudr -Bt*dudr1 + Dt*dudr2;
    case 'Steady'
        d2udrdt = 0;
end

dudr_interp = griddedInterpolant(z_u,dudr,'linear','nearest');
dudr_interp = dudr_interp(z_t);
d2udrdt_interp = griddedInterpolant(z_u,d2udrdt,'linear','nearest');
d2udrdt_interp = d2udrdt_interp(z_t);

Cc = max(sqrt((d2udrdt_interp'./dudr_interp').^2 + dudr_interp'.^2)./(L/U),1e-10)'.*a_interp.*eta_interp*Mu/SurfTens;
eta0 = (1-phi_interp).^(-1);
etainf = (1-phi_interp).^(5/3);
etar = eta_interp.*(etainf + (eta0-etainf)./(1+(6/5*Cc).^2));
etar(etar>1e12/Mu) = 1e12/Mu;

detadz = -(1./h2).*[etar(1) etar(2:2:end-2) etar(end-1)] + (1./h2).*[etar(2) etar(4:2:end) etar(end)];

etar = etar(1:2:end);

%detadz = [-A(1)*etar(1) + B(1)*etar(2) - C(1)*etar(3), ...
%          -D(2:end-1).*etar(1:end-2) - E(2:end-1).*etar(2:end-1) + C(2:end-1).*etar(3:end), ...
%          D(end)*etar(end-2) - B(end)*etar(end-1) + F(end)*etar(end)];

UP = -(1./rho_interp').*dPdz_u;
UR = (1./rho_interp').*(2*rho_interp'*g - dPdz_u*P0');
UU = 4/3*(1./rho_interp)'.*(etar'.*d2udz2 + detadz'.*dudz - etar'.*diag(4/R.^2*ones(size(u))));

switch timescheme
    case 'BDF1'
    M = eye(size([PP, PU; UP, UU])) - dt1/Re.*[PP, PU; UP, UU];
    b = ([P1 u1]' + dt1/Re.*[PR; UR]);

    case 'BDF2'
        M = eye(size([PP, PU; UP, UU])) - 1/Ft/Re.*[PP, PU; UP, UU];
        b = (-Dt/Ft*[P2 u2]' + Bt/Ft*[P1 u1]' + 1/Ft/Re.*[PR; UR]);

    case 'Steady'
        M = -[PP, PU; UP, UU];
        b = [PR; UR];
end

% BCs
% Velocity at bottom
switch BC
    case 'No normal'
        M(length(P)+1,1:end) = 0;
        M(length(P)+1,length(P)+1) = 1;
        b(length(P)+1) = 0;

    case 'No stress'
        M(length(P)+1,1:end) = 0;
        M(length(P)+1,length(P)+1:length(P)+2) = [-1, 1];
        b(length(P)+1) = 0;
end

% Pressure at surface
M(end,:) = 0;

UP_end = 0*z_p;
UP_end(end) = 1./rho_interp(end)./(z_u(end)-z_p(end));
UR_end = -dPdz_u*P0';
UR_end = UR_end(end) + (1./rho_interp(end)).*(-1./(z_u(end)-z_p(end)).*P0(end) + 2*rho_interp(end)*g);
UU_end = 0*z_u;
UU_end(end-3:end) = 1./rho_interp(end)*4/3*[-D(end).*D(end-2)*etar(end-2), -D(end)*E(end-2)*etar(end-2) + B(end)*D(end-1)*etar(end-1), ...
    C(end-2)*D(end)*etar(end-2)+B(end)*E(end-1)*etar(end-1), -B(end)*C(end-1)*etar(end-1)];

switch timescheme
    case 'BDF1'
    M(end,:) = [UP_end, UU_end] - dt1/Re.*[UP_end, UU_end];
    b(end) = (u1(end) + dt1/Re.*UR_end);

    case 'BDF2'
        M(end,:) = [UP_end, UU_end] - 1/Ft/Re.*[UP_end, UU_end];
        b(end) = -Dt/Ft*u2(end) + Bt/Ft*u1(end) + 1/Ft/Re.*UR_end;

    case 'Steady'
        M(end,:) = -[UP_end, UU_end];
        b(end) = UR_end;
end

% No strain rate at outer edge
%M(end,:) = 0;
%M(end,end-1:end) = [-1 1];
%b(end) = 0;

Vi = V;
V = M\b;
P = V(1:length(P1))';
u = V(length(P1)+1:end)';

n = n+1;

end

P = P/(L/Mu/U);
u = u*U;

%%%%%%%%%%%%% 2D Cylindrical %%%%%%%%%%%%%%%%%
function [P,u] = Cylindrical2D_pde(P1,u1,P2,u2,P0,rho,drhodt,...
            phi,a,SurfTens,eta,beta,dbetadt,R,g,z_p,z_u,dt1,dt2,BC,timescheme)

Mu = max(eta);
L = max(z_u);
Rho = max(rho);
U = max(1e-12,max(u1));
Re = U*L*Rho/Mu;

P1 = P1*(L/Mu/U);
u1 = u1/U;
P2 = P2*(L/Mu/U);
u2 = u2/U;
P0 = P0*(L/Mu/U);
rho = rho/Rho;
drhodt = drhodt/Rho*(L/U);
eta = eta/Mu;
beta = beta/(L/Mu/U);
z_p = z_p/L;
z_u = z_u/L;
dt1 = dt1/(L/U);
dt2 = dt2/(L/U);

%[~,~,~,Bt,~,Dt,~,Ft] = FDcoeff([0,dt2,dt1+dt2]);
Bt = (dt1 + dt2)/dt1/dt2;
Dt = dt1/dt2/(dt1+dt2);
Ft = (dt2+2*dt1)/dt1/(dt2+dt1);

P = P1; 
u = u1;
Vi = [P u]';
V = 0;

max_iter = 1000;
tol = 1e-7;
n = 1;

n_magmap = sqrt(length(z_p)/2);
n_magmau = sqrt(length(z_u)/2);

z_pr = z_p(1:n_magmap^2);
z_pz = z_p(n_magmap^2+1:end);%/L;
dp = (sqrt((z_pr-z_pr').^2 + (z_pz-z_pz').^2)).*(triu(ones(n_magmap^2),1) - triu(ones(n_magmap^2),1)');

z_ur = z_u(1:n_magmau^2);
z_uz = z_u(n_magmau^2+1:end);%/L;
du = (sqrt((z_ur-z_ur').^2 + (z_uz-z_uz').^2)).*(triu(ones(n_magmau^2),1) - triu(ones(n_magmau^2),1)');


% Construct first derivative in dimension 1
[h1,h2,A,B,C,D,E,F] = FDcoeff(reshape([[0; diag(dp,-1)] diag(dp) [diag(dp,1); 0]]',1,[]));

h1 = h1(2:3:end); h1(1:n_magmap:end) = h1(2:n_magmap:end);
h2 = h2(2:3:end); h2(n_magmap:n_magmap:end) = h2(n_magmap-1:n_magmap:end);

A = A(2:3:end); A(1:n_magmap:end) = A(2:n_magmap:end); A(n_magmap:n_magmap:end) = A(n_magmap-1:n_magmap:end);
B = B(2:3:end); B(1:n_magmap:end) = B(2:n_magmap:end); B(n_magmap:n_magmap:end) = B(n_magmap-1:n_magmap:end);
C = C(2:3:end); C(1:n_magmap:end) = C(2:n_magmap:end); C(n_magmap:n_magmap:end) = C(n_magmap-1:n_magmap:end);
D = D(2:3:end); D(1:n_magmap:end) = D(2:n_magmap:end); D(n_magmap:n_magmap:end) = D(n_magmap-1:n_magmap:end);
E = E(2:3:end); E(1:n_magmap:end) = E(2:n_magmap:end); E(n_magmap:n_magmap:end) = E(n_magmap-1:n_magmap:end);
F = F(2:3:end); F(1:n_magmap:end) = F(2:n_magmap:end); F(n_magmap:n_magmap:end) = F(n_magmap-1:n_magmap:end);

dPdz = diag(-D(2:end),-1) + diag(-E) + diag(C(1:end-1),1);
dPdz(1,1:3) = [-A(1), B(1), -C(1)];
dPdz(end,end-2:end) = [D(end), -B(end), F(end)];

dP0dz = [-A(1)*P0(1) + B(1)*P0(2) - C(1)*P0(3), ...
        -D(2:end-1).*P0(1:end-2) - E(2:end-1).*P0(2:end-1) + C(2:end-1).*P0(3:end), ...
        D(end).*P0(end-2) - B(end).*P0(end-1) + F(end).*P0(end)];
drhodz = [-A(1)*rho(2) + B(1)*rho(4) - C(1)*rho(6), ...
        -D(2:end-1).*rho(2:2:end-4) - E(2:end-1).*rho(4:2:end-2) + C(2:end-1).*rho(6:2:end), ...
        D(end).*rho(end-5) - B(end).*rho(end-3) + F(end).*rho(end-1)];
dbetadz = [-A(1)*beta(1) + B(1)*beta(2) - C(1)*beta(3), ...
        -D(2:end-1).*beta(1:end-2) - E(2:end-1).*beta(2:end-1) + C(2:end-1).*beta(3:end), ...
        D(end).*beta(end-2) - B(end).*beta(end-1) + F(end).*beta(end)];


h2 = z_u(2:end)-z_u(1:end-1);
dudz_P = diag(-1./[h2 h2(end)]) + diag(1./h2,1);
dudz_P = dudz_P(1:end-1,1:end);

h1 = z_u(2:end-1)-z_u(1:end-2);
h2 = z_u(3:end)-z_u(2:end-1);
w1 = (z_p - z_u(1:end-1))./[h1, h2(end)];
w2 = (z_u(2:end)-z_p)./[h1, h2(end)];
u_interp = diag([w1,w1(end)]) + diag(w2,1);
u_interp = u_interp(1:end-1,:);

[h1,h2,A,B,C,D,E,F] = FDcoeff(z_u);

dudz = diag(-D(2:end),-1) + diag(-E) + diag(C(1:end-1),1);
dudz(1,1:3) = [-A(1), B(1), -C(1)];
dudz(end,end-2:end) = [D(end), -B(end), F(end)];

h1 = z_u(2:end-1)-z_u(1:end-2);
h1 = [h1(1),h1,h1(end)];
h2 = z_u(3:end)-z_u(2:end-1);
h2 = [h2(1),h2,h2(end)];

d2udz2 = diag(2*h2(2:end)./(h1(2:end).*h2(2:end).*(h1(2:end)+h2(2:end))),-1) +...
    diag(-2*(h1+h2)./(h1.*h2.*(h1+h2))) + ...
    diag(2*h1(1:end-1)./(h1(1:end-1).*h2(1:end-1).*(h1(1:end-1)+h2(1:end-1))),1);
d2udz2(1,1:4) = [2*(3*h1(1) + 2*h2(1) + h2(3))./h1(1)./(h1(1) + h2(1))./(h1(1)+h2(1)+h2(3)), ...
                 -2*(2*h1(1) + 2*h2(1) + h2(3))./h1(1)./h2(1)./(h2(1) + h2(3)),...
                 2*(2*h1(1) + h2(1) + h2(3))./(h1(1) + h2(1))./h2(1)./h2(3),...
                 -2*(2*h1(1) + h2(1))./(h1(1) + h2(1) + h2(3))./(h2(1) + h2(3))./h2(3)];
d2udz2(end,end-3:end) = [-2*(h2(end-2) + 2*h2(end))./h1(end-2)./(h1(end-2)+h2(end-2))./(h1(end-2) + h2(end-2) + h2(end)),...
                         2*(h1(end-2) + h2(end-2) + 2*h2(end))./h1(end-2)./h2(end-2)./(h2(end-2)+h2(end)),...
                         -2*(h1(end-2) + 2*h2(end-2) + 2*h2(end))./(h1(end-2)+h2(end-2))./h2(end-2)./h2(end),...
                         2*(h1(end-2) + 2*h2(end-2) + 3*h2(end))./(h1(end-2) + h2(end-2) + h2(end))./(h2(end-2) + h2(end))./h2(end)];

h2 = [z_p(1) z_p(2:end) - z_p(1:end-1) z_u(end)-z_p(end)];
dPdz_u = diag(-1./h2(2:end),-1) + diag(1./h2);
dPdz_u(1,1:2) = dPdz_u(2,1:2);
dPdz_u = dPdz_u(:,1:end-1);

z_t = reshape([z_u; z_p, 0],1,[]);
z_t = z_t(1:end-1);

a_interp = griddedInterpolant(z_p,a,'linear','nearest');
a_interp = a_interp(z_t);
phi_interp = griddedInterpolant(z_p,phi,'linear','nearest');
phi_interp = phi_interp(z_t);

eta_interp = eta;
rho_interp = rho(1:2:end);

while (n<max_iter && norm(V-Vi)>tol) || n==1

PP = diag(-1./rho(2:2:end).*(drhodt + 1/2.*drhodz.*(u_interp*u')')) + ...
     diag(-1./beta.*(dbetadt + 1/2.*dbetadz.*(u_interp*u')')) + ...
     (-1/2*(u_interp*u').*(dPdz)) + ...
     diag(-1/2*(dudz_P*u')');
PR = (-1./rho(2:2:end).*(drhodt).*(1./beta-P0))' +...
    ((-1./beta).*dbetadt.*(-P0))';

PU = (-1./rho(2:2:end).*drhodz.*(1./beta + 1/2*P-P0))'.*u_interp + ...
    (-1./beta.*dbetadz.*(1/2*P-P0))'.*u_interp + ...
    (-1/2*(dPdz*P') + dP0dz').*u_interp + ...
    (-(1./beta + 1/2*P-P0))'.*dudz_P;

dudr = -3*u'./R;
dudr1 = -3*u1'./R;
dudr2 = -3*u2'./R;

switch timescheme
    case 'BDF1'
        d2udrdt = 0*dudr;
    case 'BDF2'
        d2udrdt = Ft*dudr -Bt*dudr1 + Dt*dudr2;
    case 'Steady'
        d2udrdt = 0;
end

dudr_interp = griddedInterpolant(z_u,dudr,'linear','nearest');
dudr_interp = dudr_interp(z_t);
d2udrdt_interp = griddedInterpolant(z_u,d2udrdt,'linear','nearest');
d2udrdt_interp = d2udrdt_interp(z_t);

Cc = max(sqrt((d2udrdt_interp'./dudr_interp').^2 + dudr_interp'.^2)./(L/U),1e-10)'.*a_interp.*eta_interp*Mu/SurfTens;
eta0 = (1-phi_interp).^(-1);
etainf = (1-phi_interp).^(5/3);
etar = eta_interp.*(etainf + (eta0-etainf)./(1+(6/5*Cc).^2));
etar(etar>1e12/Mu) = 1e12/Mu;

detadz = -(1./h2).*[etar(1) etar(2:2:end-2) etar(end-1)] + (1./h2).*[etar(2) etar(4:2:end) etar(end)];

etar = etar(1:2:end);

%detadz = [-A(1)*etar(1) + B(1)*etar(2) - C(1)*etar(3), ...
%          -D(2:end-1).*etar(1:end-2) - E(2:end-1).*etar(2:end-1) + C(2:end-1).*etar(3:end), ...
%          D(end)*etar(end-2) - B(end)*etar(end-1) + F(end)*etar(end)];

UP = -(1./rho_interp').*dPdz_u;
UR = (1./rho_interp').*(2*rho_interp'*g - dPdz_u*P0');
UU = 4/3*(1./rho_interp)'.*(etar'.*d2udz2 + detadz'.*dudz - etar'.*diag(4/R.^2*ones(size(u))));

switch timescheme
    case 'BDF1'
    M = eye(size([PP, PU; UP, UU])) - dt1/Re.*[PP, PU; UP, UU];
    b = ([P1 u1]' + dt1/Re.*[PR; UR]);

    case 'BDF2'
        M = eye(size([PP, PU; UP, UU])) - 1/Ft/Re.*[PP, PU; UP, UU];
        b = (-Dt/Ft*[P2 u2]' + Bt/Ft*[P1 u1]' + 1/Ft/Re.*[PR; UR]);

    case 'Steady'
        M = -[PP, PU; UP, UU];
        b = [PR; UR];
end

% BCs
% Velocity at bottom
switch BC
    case 'No normal'
        M(length(P)+1,1:end) = 0;
        M(length(P)+1,length(P)+1) = 1;
        b(length(P)+1) = 0;

    case 'No stress'
        M(length(P)+1,1:end) = 0;
        M(length(P)+1,length(P)+1:length(P)+2) = [-1, 1];
        b(length(P)+1) = 0;
end

% Pressure at surface
M(end,:) = 0;

UP_end = 0*z_p;
UP_end(end) = 1./rho_interp(end)./(z_u(end)-z_p(end));
UR_end = -dPdz_u*P0';
UR_end = UR_end(end) + (1./rho_interp(end)).*((-1./(z_u(end)-z_p(end))+F(end)).*P0(end) + 2*rho_interp(end)*g);
UU_end = 0*z_u;
UU_end(end-3:end) = 1./rho_interp(end)*4/3*[-D(end).^2*etar(end-2), -D(end)*E(end)*etar(end-2) + B(end)*D(end)*etar(end-1), ...
    C(end)*D(end)*etar(end-2)+B(end)*E(end)*etar(end-1), -B(end)*C(end)*etar(end-1)];

switch timescheme
    case 'BDF1'
    M(end,:) = [UP_end, UU_end] - dt1/Re.*[UP_end, UU_end];
    b(end) = (u1(end) + dt1/Re.*UR_end);

    case 'BDF2'
        M(end,:) = [UP_end, UU_end] - 1/Ft/Re.*[UP_end, UU_end];
        b(end) = -Dt/Ft*u2(end) + Bt/Ft*u1(end) + 1/Ft/Re.*UR_end;

    case 'Steady'
        M(end,:) = -[UP_end, UU_end];
        b(end) = UR_end;
end

Vi = V;
V = M\b;
P = V(1:length(P1))';
u = V(length(P1)+1:end)';

n = n+1;

end

P = P/(L/Mu/U);
u = u*U;

%%%%%%%%%%%%% Spherical %%%%%%%%%%%%%%%%%
function [P,u] = Spherical_pde(P1,u1,P2,u2,P0,rho,drhodt,...
            phi,a,SurfTens,eta,beta,dbetadt,R,g,z_p,z_u,dt1,dt2,BC,timescheme)

Mu = max(eta);
L = max(z_u);
Rho = max(rho);
U = max(1e-12,max(u1));
Re = U*L*Rho/Mu;

P1 = P1*(L/Mu/U);
u1 = u1/U;
P2 = P2*(L/Mu/U);
u2 = u2/U;
P0 = P0*(L/Mu/U);
rho = rho/Rho;
drhodt = drhodt/Rho*(L/U);
eta = eta/Mu;
beta = beta/(L/Mu/U);
z_p = z_p/L;
z_u = z_u/L;
dt1 = dt1/(L/U);
dt2 = dt2/(L/U);

%[~,~,~,Bt,~,Dt,~,Ft] = FDcoeff([0,dt2,dt1+dt2]);

Bt = (dt1 + dt2)/dt1/dt2;
Dt = dt1/dt2/(dt1+dt2);
Ft = (dt2+2*dt1)/dt1/(dt2+dt1);

P = P1;
u = u1;

Vi = [P u]';
V = 0;

max_iter = 200;
tol = 1e-7;
n = 1;

[h1,h2,A,B,C,D,E,F] = FDcoeff(z_p);
dPdz = diag(-D(2:end),-1) + diag(-E) + diag(C(1:end-1),1);
dPdz(1,1:3) = [-A(1), B(1), -C(1)];
dPdz(end,end-2:end) = [D(end), -B(end), F(end)];

dP0dz = [-A(1)*P0(1) + B(1)*P0(2) - C(1)*P0(3), ...
        -D(2:end-1).*P0(1:end-2) - E(2:end-1).*P0(2:end-1) + C(2:end-1).*P0(3:end), ...
        D(end).*P0(end-2) - B(end).*P0(end-1) + F(end).*P0(end)];
drhodz = [-A(1)*rho(2) + B(1)*rho(4) - C(1)*rho(6), ...
        -D(2:end-1).*rho(2:2:end-4) - E(2:end-1).*rho(4:2:end-2) + C(2:end-1).*rho(6:2:end), ...
        D(end).*rho(end-5) - B(end).*rho(end-3) + F(end).*rho(end-1)];
dbetadz = [-A(1)*beta(1) + B(1)*beta(2) - C(1)*beta(3), ...
        -D(2:end-1).*beta(1:end-2) - E(2:end-1).*beta(2:end-1) + C(2:end-1).*beta(3:end), ...
        D(end).*beta(end-2) - B(end).*beta(end-1) + F(end).*beta(end)];


h2 = z_u(2:end)-z_u(1:end-1);
dudz_P = diag(-1./[h2 h2(end)]) + diag(1./h2,1);
dudz_P = dudz_P(1:end-1,1:end);

h1 = z_u(2:end-1)-z_u(1:end-2);
h2 = z_u(3:end)-z_u(2:end-1);
w1 = (z_p - z_u(1:end-1))./[h1, h2(end)];
w2 = (z_u(2:end)-z_p)./[h1, h2(end)];
u_interp = diag([w1,w1(end)]) + diag(w2,1);
u_interp = u_interp(1:end-1,:);

[h1,h2,A,B,C,D,E,F] = FDcoeff(z_u);
dudz = diag(-D(2:end),-1) + diag(-E) + diag(C(1:end-1),1);
dudz(1,1:3) = [-A(1), B(1), -C(1)];
dudz(end,end-2:end) = [D(end), -B(end), F(end)];

h1 = z_u(2:end-1)-z_u(1:end-2);
h1 = [h1(1),h1,h1(end)];
h2 = z_u(3:end)-z_u(2:end-1);
h2 = [h2(1),h2,h2(end)];

z_u_fin = z_u;
z_u_fin(1) = z_u_fin(2)/1e3;

d2udz2 = diag(2*h2(2:end)./(h1(2:end).*h2(2:end).*(h1(2:end)+h2(2:end))),-1) +...
    diag(-2*(h1+h2)./(h1.*h2.*(h1+h2))) + ...
    diag(2*h1(1:end-1)./(h1(1:end-1).*h2(1:end-1).*(h1(1:end-1)+h2(1:end-1))),1);
d2udz2(1,1:4) = [2*(3*h1(1) + 2*h2(1) + h2(3))./h1(1)./(h1(1) + h2(1))./(h1(1)+h2(1)+h2(3)), ...
                 -2*(2*h1(1) + 2*h2(1) + h2(3))./h1(1)./h2(1)./(h2(1) + h2(3)),...
                 2*(2*h1(1) + h2(1) + h2(3))./(h1(1) + h2(1))./h2(1)./h2(3),...
                 -2*(2*h1(1) + h2(1))./(h1(1) + h2(1) + h2(3))./(h2(1) + h2(3))./h2(3)];
d2udz2(end,end-3:end) = [-2*(h2(end-2) + 2*h2(end))./h1(end-2)./(h1(end-2)+h2(end-2))./(h1(end-2) + h2(end-2) + h2(end)),...
                         2*(h1(end-2) + h2(end-2) + 2*h2(end))./h1(end-2)./h2(end-2)./(h2(end-2)+h2(end)),...
                         -2*(h1(end-2) + 2*h2(end-2) + 2*h2(end))./(h1(end-2)+h2(end-2))./h2(end-2)./h2(end),...
                         2*(h1(end-2) + 2*h2(end-2) + 3*h2(end))./(h1(end-2) + h2(end-2) + h2(end))./(h2(end-2) + h2(end))./h2(end)];


h2 = [z_p(1) z_p(2:end) - z_p(1:end-1) z_u(end)-z_p(end)];
dPdz_u = diag(-1./h2(2:end),-1) + diag(1./h2);
dPdz_u(1,1:2) = dPdz_u(2,1:2);
dPdz_u = dPdz_u(:,1:end-1);

z_t = reshape([z_u; z_p, 0],1,[]);
z_t = z_t(1:end-1);

a_interp = griddedInterpolant(z_p,a,'linear','nearest');
a_interp = a_interp(z_t);
phi_interp = griddedInterpolant(z_p,phi,'linear','nearest');
phi_interp = phi_interp(z_t);

eta_interp = eta;
rho_interp = rho(1:2:end);

while (n<max_iter && norm(V-Vi)>tol) || n==1

PP = diag(-1./rho(2:2:end).*(drhodt + 1/2.*drhodz.*(u_interp*u')')) + ...
     diag(-1./beta.*(dbetadt + 1/2.*dbetadz.*(u_interp*u')')) + ...
     (-1/2*(u_interp*u').*(dPdz)) + ...
     diag(-1/2*(dudz_P*u')') + ...
     (-1./z_p.*(u_interp*u')').*eye(length(z_p));

PR = (-1./rho(2:2:end).*(drhodt).*(1./beta-P0))' + ...
    (-1./beta.*(dbetadt).*(-P0))';

PU = (-1./rho(2:2:end).*drhodz.*(1./beta + 1/2*P-P0))'.*u_interp + ...
    (-1./beta.*dbetadz.*(1/2*P-P0))'.*u_interp + ...
    (-1/2*(dPdz*P') + dP0dz').*u_interp + ...
    (-(1./beta + 1/2*P-P0))'.*dudz_P + ...
    (-1./z_p.*(1./beta + 1/2*P-P0))'.*u_interp;

dudr = dudz*u'; 
dudr1 = dudz*u1';
dudr2 = dudz*u2';

switch timescheme
    case 'BDF1'
        d2udrdt = (dudr1 - dudr2)./dt1;
    case 'BDF2'
        if n == 1
            d2udrdt = (dudr1 - dudr2)./dt1;
        else
            d2udrdt = Ft.*dudr - Bt.*dudr1 + Dt.*dudr2;
        end
    case 'Steady'
        d2udrdt = 0;
end

dudr_interp = griddedInterpolant(z_u,dudr,'linear','nearest');
dudr_interp = dudr_interp(z_t);
d2udrdt_interp = griddedInterpolant(z_u,d2udrdt,'linear','nearest');
d2udrdt_interp = d2udrdt_interp(z_t);

Cc = max(sqrt((d2udrdt_interp'./dudr_interp').^2 + dudr_interp'.^2)./(L/U),1e-10)'.*a_interp.*eta_interp*Mu/SurfTens;
eta0 = (1-phi_interp).^(-1);
etainf = (1-phi_interp).^(5/3);
etar = eta_interp.*(etainf + (eta0-etainf)./(1+(6/5*Cc).^2));
etar(etar>1e12/Mu) = 1e12/Mu;

detadz = -(1./h2).*[etar(1) etar(2:2:end-2) etar(end-2)] + (1./h2).*[etar(2) etar(4:2:end) etar(end-1)];

etar = etar(1:2:end);

UP = -(1./rho_interp').*dPdz_u;
UR = 0*z_u';
UU = 4/3*(1./rho_interp)'.*(etar'.*d2udz2 + 2.*etar'./z_u_fin'.*dudz - 2.*etar'./(z_u_fin'.^2).*eye(length(u)) + ...
    detadz'.*dudz - detadz'.*(1./(z_u_fin').*eye(length(u))));

switch timescheme
    case 'BDF1'
    M = eye(size([PP, PU; UP, UU])) - dt1/Re.*[PP, PU; UP, UU];
    b = ([P1 u1]' + dt1/Re.*[PR; UR]);

    case 'BDF2'
        M = eye(size([PP, PU; UP, UU])) - 1/Ft/Re.*[PP, PU; UP, UU];
        b = (-Dt/Ft.*[P2 u2]' + Bt/Ft.*[P1 u1]' + 1/Ft/Re.*[PR; UR]);

    case 'Steady'
        M = -[PP, PU; UP, UU];
        b = [PR; UR];
end

% BCs

% Symmetry in velocity
M(length(P)+1,:) = 0;
M(length(P)+1,length(P)+1) = 1;
b(length(P)+1) = 0;

% % Pressure at surface
M(end,:) = 0;

UP_end = 0*z_p;
UP_end(end) = 1./rho_interp(end)./(z_u(end)-z_p(end));
UR_end = 0*P0';
UR_end = UR_end(end) + (1./rho_interp(end)).*(-1./(z_u(end)-z_p(end)).*P0(end));
UU_end = 0*z_u;
UU_end(end-3:end) = 1./rho_interp(end)*4/3*[(z_u(end-2)./z_u(end)).^2.*D(end).*etar(end-2).*(-D(end-2)),...
    (z_u(end-2)./z_u(end)).^2.*D(end)*etar(end-2).*(-E(end-2)-1./z_u(end-2)) + (z_u(end-1)./z_u(end)).^2.*(-B(end))*etar(end-1).*(-D(end-1)) + etar(end).*D(end),...
    (z_u(end-2)./z_u(end)).^2.*D(end)*etar(end-2).*(C(end-2)) + (z_u(end-1)./z_u(end)).^2.*(-B(end))*etar(end-1).*(-E(end-1)-1./z_u(end-1)) + etar(end).*(-B(end)),...
    (z_u(end-1)./z_u(end)).^2.*(-B(end))*etar(end-1).*C(end-1) + etar(end).*(F(end) - 1./z_u(end))];

switch timescheme
    case 'BDF1'
    M(end,:) = [UP_end, UU_end] - dt1/Re.*[UP_end, UU_end];
    b(end) = (u1(end) + dt1/Re.*UR_end);

    case 'BDF2'
        M(end,:) = [UP_end, UU_end] - 1/Ft/Re.*[UP_end, UU_end];
        b(end) = -Dt/Ft*u2(end) + Bt/Ft*u1(end) + 1/Ft/Re.*UR_end;

    case 'Steady'
        M(end,:) = -[UP_end, UU_end];
        b(end) = UR_end;
end

%M(length(P),:) = 0;
%M(length(P),length(P)) = 1;
%b(length(P)) = P0(end);

Vi = V;
V = M\b;
P = V(1:length(P))';
%if any(P<0)
%    P = P0
%end
u = V(length(P)+1:end)';

n = n+1;

end

P = P/(L/Mu/U);
u = u*U;


function [h1,h2,A,B,C,D,E,F] = FDcoeff(z)
h1 = [z(2)-z(1), z(2:end-1)-z(1:end-2), z(end-1)-z(end-2)];
h2 = [z(3)-z(2), z(3:end)-z(2:end-1), z(end)-z(end-1)];
A = (2*h1 + h2)./h1./(h1+h2);
B = (h1+h2)./h1./h2;
C = h1./(h1+h2)./h2;
D = h2./h1./(h1+h2);
E = (h1-h2)./h1./h2; 
F = (h1 + 2*h2)./h2./(h1+h2);
