function [DynFun] = getFunctions_dynamic_dimensionless(Geometry,BC)

%warning('off','MATLAB:illConditionedMatrix')
%warning('off','MATLAB:nearlySingularMatrix')

switch Geometry
    case 'Radial'
        DynFun = @(P1,u1,P2,u2,P0,rho,drhodt,phi,a,SurfTens,eta,beta,dbetadt,...
            R,g,z_p,z_u,z_t,dt1,dt2,timescheme)Spherical_pde(P1,u1,P2,u2,P0,rho,drhodt,...
            phi,a,SurfTens,eta,beta,dbetadt,R,g,z_p,z_u,z_t,dt1,dt2,BC,timescheme);

    case 'Cylindrical'
        DynFun = @(P1,u1,P2,u2,P0,rho,drhodt,phi,a,SurfTens,eta,beta,dbetadt,...
            R,g,z_p,z_u,z_t,dt1,dt2,timescheme)Cylindrical_pde(P1,u1,P2,u2,P0,rho,drhodt,...
            phi,a,SurfTens,eta,beta,dbetadt,R,g,z_p,z_u,z_t,dt1,dt2,BC,timescheme);

    case '2D Cylindrical'
        DynFun = @(P1,u1,P2,u2,P0,rho,drhodt,phi,a,SurfTens,eta,beta,dbetadt,...
            R,g,z_p,z_u,z_t,dt1,dt2,timescheme)Cylindrical2D_pde(P1,u1,P2,u2,P0,rho,drhodt,...
            phi,a,SurfTens,eta,beta,dbetadt,R,g,z_p,z_u,z_t,dt1,dt2,BC,timescheme);
    case '2D Planar'
        DynFun = @(P1,u1,P2,u2,P0,rho,drhodt,phi,a,SurfTens,eta,beta,dbetadt,...
            R,g,z_p,z_u,z_t,dt1,dt2,timescheme)Planar2D_pde(P1,u1,P2,u2,P0,rho,drhodt,...
            phi,a,SurfTens,eta,beta,dbetadt,R,g,z_p,z_u,z_t,dt1,dt2,BC,timescheme);

end

%%%%%%%%%%%%% Cylindrical %%%%%%%%%%%%%%%%%
function [P,u] = Cylindrical_pde(P1,u1,P2,u2,P0,rho,drhodt,...
            phi,a,SurfTens,eta,beta,dbetadt,R,g,z_p,z_u,z_t,dt1,dt2,BC,timescheme)

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
z_t = z_t/L;
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

%z_t = reshape([z_u; z_p, 0],1,[]);
%z_t = z_t(1:end-1);

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
UR = (1./rho_interp').*0; %*(rho_interp'*g - dPdz_u*P0');
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
UR_end = P0(end)./rho_interp(end)./(z_u(end)-z_p(end)); %-dPdz_u*P0';
%UR_end = UR_end(end) + (1./rho_interp(end)).*(1./(z_u(end)-z_p(end)).*P0(end) + 2*rho_interp(end)*g);
%UU_end = 0*z_u;
%UU_end(end-3:end) = 1./rho_interp(end)*4/3*[-D(end).*D(end-2)*etar(end-2), -D(end)*E(end-2)*etar(end-2) + B(end)*D(end-1)*etar(end-1), ...
%    C(end-2)*D(end)*etar(end-2)+B(end)*E(end-1)*etar(end-1), -B(end)*C(end-1)*etar(end-1)];
eye0 = diag(ones(length(u)-1));
UU_end = 1./rho_interp(end).*(4/3.*(D(end)*etar(end-4)*dudz(end-2,:) + ...
    -B(end)*etar(end-2)*dudz(end-1,:) + ...
    -4*etar(end)/R^2*eye0(end,:)));


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
            phi,a,SurfTens,eta,beta,dbetadt,R,g,z_p,z_u,z_t,dt1,dt2,BC,timescheme)

% Non-dimensionalization
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
z_t = z_t/L;
dt1 = dt1/(L/U);
dt2 = dt2/(L/U);

% Time stepping scheme
Bt = (dt1 + dt2)/dt1/dt2;
Dt = dt1/dt2/(dt1+dt2);
Ft = (dt2+2*dt1)/dt1/(dt2+dt1);

% Initial conditions
P = P1; 
u = u1;
Vi = [P u]';
V = 0;

% Numerical tolerance
max_iter = 100;
tol = 1e-4;
n = 1;

% Split coordinates
n_magma = floor((sqrt(length(z_t)/2)-1)/2);

z_pr = z_p(1:n_magma^2);
z_pz = z_p(n_magma^2+1:end);

z_urr = z_u(1:n_magma*(n_magma+1));
z_urz = z_u(n_magma*(n_magma+1)+1:2*n_magma*(n_magma+1));
z_uzr = z_u(2*n_magma*(n_magma+1)+1:3*n_magma*(n_magma+1));
z_uzz = z_u(3*n_magma*(n_magma+1)+1:end);

% Calculate derivatives for variables
[ddr_p,ddz_p,d2dr_p,d2dz_p] = Derivs2Diso(z_p,n_magma,n_magma);
[ddr_t,ddz_t,d2dr_t,d2dz_t] = Derivs2D(z_t,2*n_magma+1,2*n_magma+1);
[ddr_ur,ddz_ur,d2dr_ur,d2dz_ur] = Derivs2D(z_u(1:2*n_magma*(n_magma+1)),n_magma+1,n_magma);
[ddr_uz,ddz_uz,d2dr_uz,d2dz_uz] = Derivs2D(z_u(2*n_magma*(n_magma+1)+1:end),n_magma,n_magma+1);

% Calculate derivative and interpolation of velocity on staggered nodes
% Radial
I = ones(size(z_urr));
I(n_magma+1:n_magma+1:end) = 0;

h1 = z_urr(2:end)-z_urr(1:end-1);
w1 = (z_pr - z_urr(logical(I)))./h1(logical(I(1:end-1)));
w2 = (z_urr(logical(circshift(I,1))) - z_pr)./h1(logical(I(1:end-1)));

ur_interp = zeros([n_magma^2,n_magma*(n_magma+1)]);
durdr_P = zeros([n_magma^2,n_magma*(n_magma+1)]);
for i = 1:n_magma
    durdr_P((i-1)*n_magma+1:i*n_magma,(i-1)*(n_magma+1)+1:i*(n_magma+1)-1) = -diag(1./h1((i-1)*(n_magma+1)+1:i*(n_magma+1)-1));
    durdr_P((i-1)*n_magma+1:i*n_magma,(i-1)*(n_magma+1)+2:i*(n_magma+1)) = durdr_P((i-1)*n_magma+1:i*n_magma,(i-1)*(n_magma+1)+2:i*(n_magma+1)) + ...
        diag(1./h1((i-1)*(n_magma+1)+1:i*(n_magma+1)-1));
    ur_interp((i-1)*n_magma+1:i*n_magma,(i-1)*(n_magma+1)+1:i*(n_magma+1)-1) = diag(w1((i-1)*n_magma+1:i*n_magma));
    ur_interp((i-1)*n_magma+1:i*n_magma,(i-1)*(n_magma+1)+2:i*(n_magma+1)) = ur_interp((i-1)*n_magma+1:i*n_magma,(i-1)*(n_magma+1)+2:i*(n_magma+1)) + ...
        diag(w2((i-1)*n_magma+1:i*n_magma));
end

% Vertical
z_uzz_trans = reshape(reshape(z_uzz,n_magma,n_magma+1)',1,[]);
h1 = z_uzz_trans(2:end)-z_uzz_trans(1:end-1);
h1 = reshape(reshape(h1(logical(I)),n_magma,n_magma)',1,[]);

w1 = reshape(reshape(reshape(reshape(z_pz,n_magma,n_magma)',1,[]) - z_uzz_trans(logical(I)),n_magma,n_magma)',1,[])./h1;
w2 = reshape(reshape(z_uzz_trans(logical(circshift(I,1))) - reshape(reshape(z_pz,n_magma,n_magma)',1,[]),n_magma,n_magma)',1,[])./h1;

uz_interp = zeros([n_magma^2,n_magma*(n_magma+1)]);
duzdz_P = zeros([n_magma^2,n_magma*(n_magma+1)]);
uz_interp(:,1:n_magma^2) = diag(w1);
uz_interp(:,n_magma+1:end) = uz_interp(:,n_magma+1:end) + diag(w2);
duzdz_P(:,1:n_magma^2) = diag(-1./h1);
duzdz_P(:,n_magma+1:end) = duzdz_P(:,n_magma+1:end) + diag(1./h1);

% Calculate derivatives on velocity nodes
% Radial

h1 = 0*z_urr;
h1(logical(circshift(I,1))) = [z_pr(2:end) - z_pr(1:end-1), 0];
h1(1:n_magma+1:end) = z_pr(1:n_magma:end);
h1(n_magma+1:n_magma+1:end) = z_urr(n_magma+1:n_magma+1:end) - z_pr(n_magma:n_magma:end);

dPdr_ur = zeros(n_magma*(n_magma+1),n_magma^2);
for i = 1:n_magma
    dPdr_ur((i-1)*(n_magma+1)+2:i*(n_magma+1)-1,(i-1)*(n_magma)+1:i*(n_magma)-1) = diag(-1./h1((i-1)*(n_magma+1)+2:i*(n_magma+1)-1));
    dPdr_ur((i-1)*(n_magma+1)+2:i*(n_magma+1)-1,(i-1)*(n_magma)+2:i*(n_magma)) = dPdr_ur((i-1)*(n_magma+1)+2:i*(n_magma+1)-1,(i-1)*(n_magma)+2:i*(n_magma)) + ...
        diag(1./h1((i-1)*(n_magma+1)+2:i*(n_magma+1)-1));
end
dPdr_ur(1:n_magma+1:end,:) = dPdr_ur(2:n_magma+1:end,:);
dPdr_ur(n_magma+1:n_magma+1:end,:) = dPdr_ur(n_magma:n_magma+1:end,:);

% Vertical

z_pz_trans = reshape(reshape(z_pz,n_magma,[])',1,[]);
h1 = reshape(reshape([z_pz(1) z_pz_trans(2:end) - z_pz_trans(1:end-1)],n_magma,[])',1,[]);
h1(1:n_magma) = z_pz(1:n_magma);
h1 = [h1 z_uzz(end-n_magma+1:end) - z_pz(end-n_magma+1:end)];

dPdz_uz = zeros(n_magma*(n_magma+1),n_magma^2);
dPdz_uz(n_magma+1:end-n_magma,1:end-n_magma) = diag(-1./h1(n_magma+1:end-n_magma)); 
dPdz_uz(n_magma+1:end-n_magma,n_magma+1:end) = dPdz_uz(n_magma+1:end-n_magma,n_magma+1:end) + diag(1./h1(n_magma+1:end-n_magma));

dPdz_uz(1:n_magma,:) = dPdz_uz(n_magma+1:2*n_magma,:);
dPdz_uz(end-n_magma+1:end,:) = dPdz_uz(end-2*n_magma+1:end-n_magma,:);

% Cross-derivatives in velocity
% Radial 

x1 = z_uzr(1:end-n_magma); y1 = z_uzz(1:end-n_magma);
x2 = circshift(z_uzr(1:end-n_magma),-1); y2 = z_uzz(2:end-n_magma+1);
x3 = circshift(z_uzr(n_magma+1:end),-1); y3 = circshift(z_uzz(n_magma+1:end),-1);
x4 = z_uzr(n_magma+1:end); y4 = z_uzz(n_magma+1:end);

x0 = z_urr(logical(circshift(I,1))); y0 = z_urz(logical(circshift(I,1)));

a_zr = (x1 - x2 + x3 - x4).*(-y1 + y4) - (-x1 + x4).*(y1 - y2 + y3 - y4);
b_zr = ((x1 - x2 + x3 - x4).*y1 - x1.*(y1 - y2 + y3 - y4) +...
    (-x1 + x2).*(-y1 + y4) - (-x1 + x4).*(-y1 + y2) + ...
    x0.*(y1 - y2 + y3 - y4) - y0.*(x1 - x2 + x3 - x4));
c_zr = (-x1 + x2).*y1 - x1.*(-y1 + y2) + x0.*(-y1 + y2) - y0.*(-x1 + x2);

m = (-b_zr + sqrt(b_zr.^2 - 4*a_zr.*c_zr))./(2*a_zr);
m(abs(a_zr)<1e-16) = -c_zr(abs(a_zr)<1e-16)./b_zr(abs(a_zr)<1e-16);
l = (x0 - x1 - (-x1 + x4).*m)./((-x1 + x2) + (x1 - x2 + x3 - x4).*m);

w1 = (z_urz(1:n_magma+1:end) - z_uzz(1:n_magma:end-n_magma-1))./(z_uzz(n_magma+1:n_magma:end) - z_uzz(1:n_magma:end-n_magma-1));
w2 = (z_uzz(n_magma+1:n_magma:end) - z_urz(1:n_magma+1:end))./(z_uzz(n_magma+1:n_magma:end) - z_uzz(1:n_magma:end-n_magma-1));
w3 = (z_urz(n_magma+1:n_magma+1:end) - z_uzz(n_magma:n_magma:end-n_magma))./(z_uzz(2*n_magma:n_magma:end) - z_uzz(n_magma:n_magma:end-n_magma));
w4 = (z_uzz(2*n_magma:n_magma:end) - z_urz(n_magma+1:n_magma+1:end))./(z_uzz(2*n_magma:n_magma:end) - z_uzz(n_magma:n_magma:end-n_magma));

uz_I = eye(length(z_uzr));
uz_interp_ur = zeros(n_magma*(n_magma+1));
uz_interp_ur(logical(circshift(I,1)),:) = ((1-l).*(1-m))'.*uz_I(1:end-n_magma,:) + ...
    ((1-l).*m)'.*uz_I(2:end-n_magma+1,:) + (l.*(1-m))'.*[uz_I(n_magma+2:end,:); uz_I(end,:)] + ...
    (l.*m)'.*uz_I(n_magma+1:end,:);

duzdr_ur = zeros(n_magma*(n_magma+1));
duzdr_ur(logical(circshift(I,1)),:) = ((1-l).*(1-m))'.*ddr_uz(1:end-n_magma,:) + ...
    ((1-l).*m)'.*ddr_uz(2:end-n_magma+1,:) + (l.*(1-m))'.*[ddr_uz(n_magma+2:end,:); ddr_uz(end,:)] + ...
    (l.*m)'.*ddr_uz(n_magma+1:end,:);

duzdz_ur = zeros(n_magma*(n_magma+1));
duzdz_ur(logical(circshift(I,1)),:) = ((1-l).*(1-m))'.*ddz_uz(1:end-n_magma,:) + ...
    ((1-l).*m)'.*ddz_uz(2:end-n_magma+1,:) + (l.*(1-m))'.*[ddz_uz(n_magma+2:end,:); ddz_uz(end,:)] + ...
    (l.*m)'.*ddz_uz(n_magma+1:end,:);

d2uzdz2_ur = zeros(n_magma*(n_magma+1));
d2uzdz2_ur(logical(circshift(I,1)),:) = ((1-l).*(1-m))'.*d2dz_uz(1:end-n_magma,:) + ...
    ((1-l).*m)'.*d2dz_uz(2:end-n_magma+1,:) + (l.*(1-m))'.*[d2dz_uz(n_magma+2:end,:); d2dz_uz(end,:)] + ...
    (l.*m)'.*d2dz_uz(n_magma+1:end,:);

uz_interp_ur(1:n_magma+1:end,:) = w1'.*uz_I(1:n_magma:end-n_magma-1,:) + w2'.*uz_I(n_magma+1:n_magma:end,:);
duzdr_ur(1:n_magma+1:end,:) = w1'.*ddr_uz(1:n_magma:end-n_magma-1,:) + w2'.*ddr_uz(n_magma+1:n_magma:end,:);
duzdr_ur(1:n_magma+1:end,:) = 0*duzdr_ur(1:n_magma+1:end,:);
duzdz_ur(1:n_magma+1:end,:) = w1'.*ddz_uz(1:n_magma:end-n_magma-1,:) + w2'.*ddz_uz(n_magma+1:n_magma:end,:);
d2uzdz2_ur(1:n_magma+1:end,:) = w1'.*d2dz_uz(1:n_magma:end-n_magma-1,:) + w2'.*d2dz_uz(n_magma+1:n_magma:end,:);

uz_interp_ur(n_magma+1:n_magma+1:end,:) = w3'.*uz_I(n_magma:n_magma:end-n_magma,:) + w4'.*uz_I(2*n_magma:n_magma:end,:);
duzdr_ur(n_magma+1:n_magma+1:end,:) = w3'.*ddr_uz(n_magma:n_magma:end-n_magma,:) + w4'.*ddr_uz(2*n_magma:n_magma:end,:);
duzdz_ur(n_magma+1:n_magma+1:end,:) = w3'.*ddz_uz(n_magma:n_magma:end-n_magma,:) + w4'.*ddz_uz(2*n_magma:n_magma:end,:);
d2uzdz2_ur(n_magma+1:n_magma+1:end,:) = w3'.*d2dz_uz(n_magma:n_magma:end-n_magma,:) + w4'.*d2dz_uz(2*n_magma:n_magma:end,:);
switch BC{2}
    case 'Confined'
        switch BC{3}
            case 'No slip'
                uz_interp_ur(n_magma+1:n_magma+1:end,:) = (w3'.*uz_I(n_magma:n_magma:end-n_magma,:) + w4'.*uz_I(2*n_magma:n_magma:end,:))/2;
                duzdr_ur(n_magma+1:n_magma+1:end,:) = (w3'.*ddr_uz(n_magma:n_magma:end-n_magma,:) + w4'.*ddr_uz(2*n_magma:n_magma:end,:))/2;
                duzdz_ur(n_magma+1:n_magma+1:end,:) = 0*duzdz_ur(n_magma+1:n_magma+1:end,:);
                d2uzdz2_ur(n_magma+1:n_magma+1:end,:) = 0*d2uzdz2_ur(n_magma+1:n_magma+1:end,:);
        end
end

d2uzdrdz_ur = ddz_ur*duzdr_ur;

% Vertical
x1 = z_urr(logical(I)); y1 = z_urz(logical(I));
x4 = x1(n_magma+1:end); y4 = y1(n_magma+1:end); 
x1 = x1(1:end-n_magma); y1 = y1(1:end-n_magma);

x2 = z_urr(logical(circshift(I,1))); y2 = z_urz(logical(circshift(I,1))); 
x3 = x2(n_magma+1:end); y3 = y2(n_magma+1:end); 
x2 = x2(1:end-n_magma); y2 = y2(1:end-n_magma);

x0 = z_uzr(n_magma+1:end-n_magma); y0 = z_uzz(n_magma+1:end-n_magma);

a_zz = (x1 - x2 + x3 - x4).*(-y1 + y4) - (-x1 + x4).*(y1 - y2 + y3 - y4);
b_zz = ((x1 - x2 + x3 - x4).*y1 - x1.*(y1 - y2 + y3 - y4) +...
    (-x1 + x2).*(-y1 + y4) - (-x1 + x4).*(-y1 + y2) + ...
    x0.*(y1 - y2 + y3 - y4) - y0.*(x1 - x2 + x3 - x4));
c_zz = (-x1 + x2).*y1 - x1.*(-y1 + y2) + x0.*(-y1 + y2) - y0.*(-x1 + x2);

m = (-b_zz + sqrt(b_zz.^2 - 4*a_zz.*c_zz))./(2*a_zz);
m(abs(a_zz)<1e-16) = -c_zz(abs(a_zz)<1e-16)./b_zz(abs(a_zz)<1e-16);
l = (x0 - x1 - (-x1 + x4).*m)./((-x1 + x2) + (x1 - x2 + x3 - x4).*m);
ur_I = eye(length(z_urr));

w1 = (z_uzr(1:n_magma) - z_urr(1:n_magma))./(z_urr(2:n_magma+1) - z_urr(1:n_magma));
w2 = (z_urr(2:n_magma+1) - z_uzr(1:n_magma))./(z_urr(2:n_magma+1) - z_urr(1:n_magma));
w3 = (z_uzr(end-n_magma+1:end) - z_urr(end-n_magma:end-1))./(z_urr(end-n_magma+1:end) - z_urr(end-n_magma:end-1));
w4 = (z_urr(end-n_magma+1:end) - z_uzr(end-n_magma+1:end))./(z_urr(end-n_magma+1:end) - z_urr(end-n_magma:end-1));

ur_interp_uz = zeros(n_magma*(n_magma+1));
durdr_uz = zeros(n_magma*(n_magma+1));
durdz_uz = zeros(n_magma*(n_magma+1));
for i = 1:n_magma-1
ur_interp_uz(i*n_magma+1:(i+1)*n_magma,:) = ((1-l((i-1)*n_magma+1:i*n_magma)).*(1-m((i-1)*n_magma+1:i*n_magma)))'.*ur_I((i-1)*(n_magma+1)+1:i*(n_magma+1)-1,:) + ...
    ((1-l((i-1)*n_magma+1:i*n_magma)).*m((i-1)*n_magma+1:i*n_magma))'.*ur_I((i-1)*(n_magma+1)+2:i*(n_magma+1),:) + ...
    (l((i-1)*n_magma+1:i*n_magma).*(1-m((i-1)*n_magma+1:i*n_magma)))'.*ur_I(i*(n_magma+1)+2:(i+1)*(n_magma+1),:) + ...
    (l((i-1)*n_magma+1:i*n_magma).*m((i-1)*n_magma+1:i*n_magma))'.*ur_I(i*(n_magma+1)+1:(i+1)*(n_magma+1)-1,:);

durdr_uz(i*n_magma+1:(i+1)*n_magma,:) = ((1-l((i-1)*n_magma+1:i*n_magma)).*(1-m((i-1)*n_magma+1:i*n_magma)))'.*ddr_ur((i-1)*(n_magma+1)+1:i*(n_magma+1)-1,:) + ...
    ((1-l((i-1)*n_magma+1:i*n_magma)).*m((i-1)*n_magma+1:i*n_magma))'.*ddr_ur((i-1)*(n_magma+1)+2:i*(n_magma+1),:) + ...
    (l((i-1)*n_magma+1:i*n_magma).*(1-m((i-1)*n_magma+1:i*n_magma)))'.*ddr_ur(i*(n_magma+1)+2:(i+1)*(n_magma+1),:) + ...
    (l((i-1)*n_magma+1:i*n_magma).*m((i-1)*n_magma+1:i*n_magma))'.*ddr_ur(i*(n_magma+1)+1:(i+1)*(n_magma+1)-1,:);

durdz_uz(i*n_magma+1:(i+1)*n_magma,:) = ((1-l((i-1)*n_magma+1:i*n_magma)).*(1-m((i-1)*n_magma+1:i*n_magma)))'.*ddz_ur((i-1)*(n_magma+1)+1:i*(n_magma+1)-1,:) + ...
    ((1-l((i-1)*n_magma+1:i*n_magma)).*m((i-1)*n_magma+1:i*n_magma))'.*ddz_ur((i-1)*(n_magma+1)+2:i*(n_magma+1),:) + ...
    (l((i-1)*n_magma+1:i*n_magma).*(1-m((i-1)*n_magma+1:i*n_magma)))'.*ddz_ur(i*(n_magma+1)+2:(i+1)*(n_magma+1),:) + ...
    (l((i-1)*n_magma+1:i*n_magma).*m((i-1)*n_magma+1:i*n_magma))'.*ddz_ur(i*(n_magma+1)+1:(i+1)*(n_magma+1)-1,:);
end

% Bottom boundary
switch BC{1}
    case 'No stress'
        ur_interp_uz(1:n_magma,:) = w1'.*ur_I(1:n_magma,:) + w2'.*ur_I(2:n_magma+1,:);
        durdr_uz(1:n_magma,:) = w1'.*ddr_ur(1:n_magma,:) + w2'.*ddr_ur(2:n_magma+1,:);
        durdz_uz(1:n_magma,:) = w1'.*ddz_ur(1:n_magma,:) + w2'.*ddz_ur(2:n_magma+1,:);
    case 'No normal'
        ur_interp_uz(1:n_magma,:) = 0;
        durdr_uz(1:n_magma,:) = 0;
        durdz_uz(1:n_magma,:) = (w1'.*ddz_ur(1:n_magma,:) + w2'.*ddz_ur(2:n_magma+1,:))/2;
end
% Top boundary
switch BC{4}
    case 'Load'
        ur_interp_uz(end-n_magma+1:end,:) = w3'.*ur_I(end-n_magma:end-1,:) + w4'.*ur_I(end-n_magma+1:end,:);
        durdr_uz(end-n_magma+1:end,:) = w3'.*ddr_ur(end-n_magma:end-1,:) + w4'.*ddr_ur(end-n_magma+1:end,:);
        durdz_uz(end-n_magma+1:end,:) = w3'.*ddz_ur(end-n_magma:end-1,:) + w4'.*ddz_ur(end-n_magma+1:end,:);
    case 'Dirichlet'
        ur_interp_uz(end-n_magma+1:end,:) = 0;
        durdr_uz(end-n_magma+1:end,:) = 0;
        durdz_uz(end-n_magma+1:end,:) = (w3'.*ddz_ur(end-n_magma:end-1,:) + w4'.*ddz_ur(end-n_magma+1:end,:))/2;
end

d2urdrdz_uz = ddr_uz*durdz_uz;

% Calculate derivatives and interpolations of constants
dP0dz = (ddz_p*P0')';
drhodr = (ddr_t*rho')';
drhodz = (ddz_t*rho')';
dbetadr = (ddz_p*beta')';
dbetadz = (ddz_p*beta')';

a_interp = scatteredInterpolant(z_pr',z_pz',a','linear','nearest');
a_interp = a_interp(z_t(1:(2*n_magma+1)^2)', z_t((2*n_magma+1)^2+1:end)')';
phi_interp = scatteredInterpolant(z_pr',z_pz',phi','linear','nearest');
phi_interp = phi_interp(z_t(1:(2*n_magma+1)^2)', z_t((2*n_magma+1)^2+1:end)')';

z_t_r = reshape(z_t(1:(2*n_magma+1)^2),2*n_magma+1,2*n_magma+1);
z_t_z = reshape(z_t((2*n_magma+1)^2+1:end),2*n_magma+1,2*n_magma+1);
ind_p = 0*z_t_r;
ind_p(2:2:end,2:2:end) = 1; ind_p = logical(repmat(reshape(ind_p,1,[]),1,2));
ind_ur = 0*z_t_r;
ind_ur(1:2:end,2:2:end) = 1; ind_ur = logical(repmat(reshape(ind_ur,1,[]),1,2));
ind_uz = 0*z_t_z;
ind_uz(2:2:end,1:2:end) = 1; ind_uz = logical(repmat(reshape(ind_uz,1,[]),1,2));

eta_interp = eta;
rho_interpr = rho(ind_ur(1:(2*n_magma+1)^2));
rho_interpz = rho(ind_uz(1:(2*n_magma+1)^2));

while (n<max_iter && norm(V-Vi)>tol) || n==1

PP = diag(-1./rho(ind_p(1:(2*n_magma+1)^2)).*(drhodt + ...
     1/2.*drhodr(ind_p(1:(2*n_magma+1)^2)).*(ur_interp*u(1:n_magma*(n_magma+1))')' + ...
     1/2.*drhodz(ind_p(1:(2*n_magma+1)^2)).*(uz_interp*u(n_magma*(n_magma+1)+1:end)')')) + ...
     diag(-1./beta.*(dbetadt + 1/2.*dbetadr.*(ur_interp*u(1:n_magma*(n_magma+1))')' + ...
     1/2.*dbetadz.*(uz_interp*u(n_magma*(n_magma+1)+1:end)')')) + ...
     (-1/2*(ur_interp*u(1:n_magma*(n_magma+1))').*ddr_p) + ...
     (-1/2*(uz_interp*u(n_magma*(n_magma+1)+1:end)').*(ddz_p)) + ...
     diag(-1/2./(z_pr).*(ur_interp*u(1:n_magma*(n_magma+1))')' + ...
     -1/2*(durdr_P*u(1:n_magma*(n_magma+1))')' + ...
     -1/2*(duzdz_P*u(n_magma*(n_magma+1)+1:end)')');
PR = (-1./rho(ind_p(1:(2*n_magma+1)^2)).*(drhodt).*(1./beta-P0))' +...
    ((-1./beta).*dbetadt.*(-P0))';
PUr = (-1./rho(ind_p(1:(2*n_magma+1)^2)).*drhodr(ind_p(1:(2*n_magma+1)^2)).*(1./beta + 1/2*P-P0))'.*ur_interp + ...
    (-1./beta.*dbetadr.*(1/2*P-P0))'.*ur_interp + ...
    (-1/2*ddr_p*P').*ur_interp + ...
    (-(1./beta + 1/2*P-P0))'.*(1./z_pr'.*ur_interp + durdr_P);
PUz = (-1./rho(ind_p(1:(2*n_magma+1)^2)).*drhodz(ind_p(1:(2*n_magma+1)^2)).*(1./beta + 1/2*P-P0))'.*uz_interp + ...
    (-1./beta.*dbetadz.*(1/2*P-P0))'.*uz_interp + ...
    (-1/2*(ddz_p*P') + dP0dz').*uz_interp + ...
    (-(1./beta + 1/2*P-P0))'.*duzdz_P;

% Calculate effective viscosity

% rinv = (1./z_urr);
% rinv(~isfinite(1./z_urr)) = 0;
% 
% gamma_ur = ddr_ur*u(1:n_magma*(n_magma+1))' + 1/2*ddz_ur*u(1:n_magma*(n_magma+1))' +...
%            (rinv.*u(1:n_magma*(n_magma+1)))' + 1/2*duzdr_ur*u(n_magma*(n_magma+1)+1:end)' + ...
%            duzdz_ur*u(n_magma*(n_magma+1)+1:end)';
% gamma_ur1 = ddr_ur*u1(1:n_magma*(n_magma+1))' + 1/2*ddz_ur*u1(1:n_magma*(n_magma+1))' +...
%            (rinv.*u1(1:n_magma*(n_magma+1)))' + 1/2*duzdr_ur*u1(n_magma*(n_magma+1)+1:end)' + ...
%            duzdz_ur*u1(n_magma*(n_magma+1)+1:end)';
% gamma_ur2 = ddr_ur*u2(1:n_magma*(n_magma+1))' + 1/2*ddz_ur*u2(1:n_magma*(n_magma+1))' +...
%            (rinv.*u2(1:n_magma*(n_magma+1)))' + 1/2*duzdr_ur*u2(n_magma*(n_magma+1)+1:end)' + ...
%            duzdz_ur*u2(n_magma*(n_magma+1)+1:end)';
% 
% gamma_uz = durdr_uz*u(1:n_magma*(n_magma+1))' + 1/2*durdz_uz*u(1:n_magma*(n_magma+1))' +...
%            (1./z_uzr)'.*(ur_interp_uz*u(1:n_magma*(n_magma+1))') + 1/2*ddr_uz*u(n_magma*(n_magma+1)+1:end)' + ...
%            ddz_uz*u(n_magma*(n_magma+1)+1:end)';
% gamma_uz1 = durdr_uz*u1(1:n_magma*(n_magma+1))' + 1/2*durdz_uz*u1(1:n_magma*(n_magma+1))' +...
%            (1./z_uzr)'.*(ur_interp_uz*u1(1:n_magma*(n_magma+1))') + 1/2*ddr_uz*u1(n_magma*(n_magma+1)+1:end)' + ...
%            ddz_uz*u1(n_magma*(n_magma+1)+1:end)';
% gamma_uz2 = durdr_uz*u2(1:n_magma*(n_magma+1))' + 1/2*durdz_uz*u2(1:n_magma*(n_magma+1))' +...
%            (1./z_uzr)'.*(ur_interp_uz*u2(1:n_magma*(n_magma+1))') + 1/2*ddr_uz*u2(n_magma*(n_magma+1)+1:end)' + ...
%            ddz_uz*u2(n_magma*(n_magma+1)+1:end)';

gamma_ur = ddz_ur*u(1:n_magma*(n_magma+1))' + uz_interp_ur*(ddr_uz*u(n_magma*(n_magma+1)+1:end)');
gamma_uz = ur_interp_uz*(ddz_ur*u(1:n_magma*(n_magma+1))') + ddr_uz*u(n_magma*(n_magma+1)+1:end)';

gamma_ur1 = ddz_ur*u1(1:n_magma*(n_magma+1))' + uz_interp_ur*(ddr_uz*u1(n_magma*(n_magma+1)+1:end)');
gamma_uz1 = ur_interp_uz*(ddz_ur*u1(1:n_magma*(n_magma+1))') + ddr_uz*u1(n_magma*(n_magma+1)+1:end)';

gamma_ur2 = ddz_ur*u2(1:n_magma*(n_magma+1))' + uz_interp_ur*(ddr_uz*u2(n_magma*(n_magma+1)+1:end)');
gamma_uz2 = ur_interp_uz*(ddz_ur*u2(1:n_magma*(n_magma+1))') + ddr_uz*u2(n_magma*(n_magma+1)+1:end)';

switch timescheme
    case 'BDF1'
        d2urdrdt = 0*gamma_ur;
        d2uzdrdt = 0*gamma_uz;
    case 'BDF2'
        d2urdrdt = Ft*gamma_ur -Bt*gamma_ur1 + Dt*gamma_ur2;
        d2uzdrdt = Ft*gamma_uz -Bt*gamma_uz1 + Dt*gamma_uz2;
    case 'Steady'
        d2urdrdt = 0*gamma_ur;
        d2uzdrdt = 0*gamma_uz;
end

gamma_interp = scatteredInterpolant([z_urr,z_uzr]',[z_urz,z_uzz]',[gamma_ur;gamma_uz],'linear','nearest');
gamma_interp = gamma_interp(z_t(1:(2*n_magma+1)^2),z_t((2*n_magma+1)^2+1:end));
d2udrdt_interp = scatteredInterpolant([z_urr,z_uzr]',[z_urz,z_uzz]',[d2urdrdt;d2uzdrdt],'linear','nearest');
d2udrdt_interp = d2udrdt_interp(z_t(1:(2*n_magma+1)^2),z_t((2*n_magma+1)^2+1:end));

Cc = max(sqrt((d2udrdt_interp'./gamma_interp').^2 + gamma_interp'.^2)./(L/U),1e-10)'.*a_interp.*eta_interp*Mu/SurfTens;
eta0 = (1-phi_interp).^(-1);
etainf = (1-phi_interp).^(5/3);
etar = eta_interp.*(etainf + (eta0-etainf)./(1+(6/5*Cc).^2));
etar(etar>1e12/Mu) = 1e12/Mu;

detadr = (ddr_t*etar')';
detadz = (ddz_t*etar')';

UrP = -(1./rho_interpr').*dPdr_ur;
UrR = 0*z_urr';
UrUr = 4/3*(1./rho_interpr)'.*(etar(ind_ur(1:(2*n_magma+1)^2))'.*(d2dr_ur + min(1./z_urr',1e4).*ddr_ur - min(1./z_urr',1e4).^2.*eye(length(z_urr))) + ...
    detadr(ind_ur(1:(2*n_magma+1)^2))'.*(ddr_ur - min(1./(2*z_urr)',1e4).*eye(length(z_urr)))) + ...
    (1./rho_interpr)'.*(etar(ind_ur(1:(2*n_magma+1)^2))'.*d2dz_ur) + ...
    (1./rho_interpr)'.*detadz(ind_ur(1:(2*n_magma+1)^2))'.*(ddz_ur);
UrUz = (1./rho_interpr)'.*(etar(ind_ur(1:(2*n_magma+1)^2))'.*(1/3*d2uzdrdz_ur) + ...
    detadr(ind_ur(1:(2*n_magma+1)^2))'.*(-2/3*duzdz_ur) + ...
    detadz(ind_ur(1:(2*n_magma+1)^2))'.*(duzdr_ur));

UzP = -(1./rho_interpz').*dPdz_uz;
UzR = (1./rho_interpz').*0; %(2*rho_interpz'*g - dPdz_uz*P0');
UzUr = (1./rho_interpz)'.*(etar(ind_uz(1:(2*n_magma+1)^2))'.*(1/3*d2urdrdz_uz + 1./(3*z_uzr)'.*durdz_uz) + ...
        detadr(ind_uz(1:(2*n_magma+1)^2))'.*durdz_uz + ...
        detadz(ind_uz(1:(2*n_magma+1)^2))'.*(-2/3*durdr_uz - 2./(3*z_uzr)'.*ur_interp_uz));
UzUz = (1./rho_interpz)'.*(etar(ind_uz(1:(2*n_magma+1)^2))'.*(4/3*d2dz_uz + ...
        1/3*1./z_uzr'.*ddr_uz + d2dr_uz) + ...
        detadr(ind_uz(1:(2*n_magma+1)^2))'.*ddr_uz + ...
        detadz(ind_uz(1:(2*n_magma+1)^2))'.*(4/3*ddz_uz));

switch timescheme
    case 'BDF1'
    M = eye(size([PP, PUr, PUz; UrP, UrUr, UrUz; UzP, UzUr, UzUz])) - dt1/Re.*[PP, PUr, PUz; UrP, UrUr, UrUz; UzP, UzUr, UzUz];
    b = ([P1 u1]' + dt1/Re.*[PR; UrR; UzR]);

    case 'BDF2'
        M = eye(size([PP, PUr, PUz; UrP, UrUr, UrUz; UzP, UzUr, UzUz])) - 1/Ft/Re.*[PP, PUr, PUz; UrP, UrUr, UrUz; UzP, UzUr, UzUz];
        b = (-Dt/Ft*[P2 u2]' + Bt/Ft*[P1 u1]' + 1/Ft/Re.*[PR; UrR; UzR]);

    case 'Steady'
        M = -[PP, PUr, PUz; UrP, UrUr, UrUz; UzP, UzUr, UzUz];
        b = [PR; UrR; UzR];
end

% BCs
% Velocity at bottom
switch BC{1}
    case 'No normal'
        M(length(P)+1:length(P)+n_magma+1,:) = 0;
        M(length(P)+1:length(P)+n_magma+1,length(P)+1:length(P)+n_magma+1) = eye(n_magma+1);
        b(length(P)+1:length(P)+n_magma+1) = 0;

        M(length(P)+length(z_urr)+1:length(P)+length(z_urr)+n_magma,:) = 0;
        M(length(P)+length(z_urr)+1:length(P)+length(z_urr)+n_magma,...
            length(P)+length(z_urr)+1:length(P)+length(z_urr)+n_magma) = eye(n_magma);
        b(length(P)+length(z_urr)+1:length(P)+length(z_urr)+n_magma) = 0;

    case 'No stress'
        %M(length(P)+1:length(P)+n_magma+1,:) = 0;
        %M(length(P)+1:length(P)+n_magma+1,length(P)+1:length(P)+n_magma+1) = eye(n_magma+1);
        %b(length(P)+1:length(P)+n_magma+1) = 0;

        M(length(P)+length(z_urr)+1:length(P)+length(z_urr)+n_magma,:) = 0;
        M(length(P)+length(z_urr)+1:length(P)+length(z_urr)+n_magma,...
            length(P)+length(z_urr)+1:length(P)+length(z_urr)+n_magma) = -eye(n_magma);
        M(length(P)+length(z_urr)+1:length(P)+length(z_urr)+n_magma,...
            length(P)+length(z_urr)+n_magma+1:length(P)+length(z_urr)+2*n_magma) = eye(n_magma);
        b(length(P)+length(z_urr)+1:length(P)+length(z_urr)+n_magma) = 0;
end

% Pressure at top surface
M(end-n_magma+1:end,:) = 0;

dr = [z_pr(end-n_magma+2) - z_pr(end-n_magma+1),...
    z_pr(end-n_magma+3:end) - z_pr(end-n_magma+1:end-2), ...
    z_pr(end) - z_pr(end-1)];
dz = [z_pz(end-n_magma+2) - z_pz(end-n_magma+1),...
    z_pz(end-n_magma+3:end) - z_pz(end-n_magma+1:end-2), ...
    z_pz(end) - z_pz(end-1)];
thetas = atan2(dz,dr);

UzP_top = zeros(n_magma,n_magma.^2);
UzP_top(:,end-n_magma+1:end) = -1./rho_interpz(end-n_magma+1:end)'./(z_uzz(end-n_magma+1:end)-z_p(end-n_magma+1:end))'.*eye(n_magma);
UzR_top = 1./rho_interpz(end-n_magma+1:end)'./(z_uzz(end-n_magma+1:end)-z_p(end-n_magma+1:end))'.*P0(end-n_magma+1:end)' - ...
    1./rho_interpz(end-n_magma+1:end)'.*(ddz_p(end-n_magma+1:end,:)*P0');

ind_top = find(ind_uz(1:(2*n_magma+1)^2));

UzUr_top = 1./rho_interpz(end-n_magma+1:end)'*2/3.*(diag(ddz_uz(end-n_magma+1:end,end-3*n_magma+1:end-n_magma)).*(etar(ind_top(end-3*n_magma+1:end-2*n_magma))'.*(-durdr_uz(end-3*n_magma+1:end-2*n_magma,:) - ...
    1./max(z_uzr(end-3*n_magma+1:end-2*n_magma),1e-6)'.*ur_interp_uz(end-3*n_magma+1:end-2*n_magma,:))) + ...
    diag(ddz_uz(end-n_magma+1:end,end-2*n_magma+1:end-n_magma)).*(etar(ind_top(end-2*n_magma+1:end-n_magma))'.*(-durdr_uz(end-2*n_magma+1:end-n_magma,:)) - ...
    1./max(z_uzr(end-2*n_magma+1:end-n_magma),1e-6)'.*ur_interp_uz(end-2*n_magma+1:end-n_magma,:))) + ...
    1./rho_interpz(end-n_magma+1:end)'.*diag(ddz_uz(end-n_magma+1:end,end-n_magma+1:end)).*(etar(ind_top(end-n_magma+1:end))'.*tan(thetas)'.*durdz_uz(end-n_magma+1:end,:)) + ...
    1./rho_interpz(end-n_magma+1:end)'.*etar(ind_top(end-n_magma+1:end))'.*(1./z_uzr(end-n_magma+1:end)'.*durdz_uz(end-n_magma+1:end,:) + d2urdrdz_uz(end-n_magma+1:end,:)) + ...
    1./rho_interpz(end-n_magma+1:end)'.*detadr(ind_top(end-n_magma+1:end))'.*durdz_uz(end-n_magma+1:end,:);

UzUz_top = 1./rho_interpz(end-n_magma+1:end)'*4/3.*(diag(ddz_uz(end-n_magma+1:end,end-3*n_magma+1:end-2*n_magma)).*(etar(ind_top(end-3*n_magma+1:end-2*n_magma))'.*(ddz_uz(end-3*n_magma+1:end-2*n_magma,:))) + ...
    diag(ddz_uz(end-n_magma+1:end,end-2*n_magma+1:end)).*(etar(ind_top(end-2*n_magma+1:end-n_magma))'.*(ddz_uz(end-2*n_magma+1:end-n_magma,:)))) + ...
    1./rho_interpz(end-n_magma+1:end)'.*diag(ddz_uz(end-n_magma+1:end,end-n_magma+1:end)).*(etar(ind_top(end-n_magma+1:end))'.*tan(thetas)'.*ddr_uz(end-n_magma+1:end,:)) + ...
    1./rho_interpz(end-n_magma+1:end)'.*etar(ind_top(end-n_magma+1:end))'.*(1./z_uzr(end-n_magma+1:end)'.*ddr_uz(end-n_magma+1:end,:) + d2dr_uz(end-n_magma+1:end,:)) + ...
    1./rho_interpz(end-n_magma+1:end)'.*detadr(ind_top(end-n_magma+1:end))'.*ddr_uz(end-n_magma+1:end,:);

switch timescheme
    case 'BDF1'
        M(end-n_magma+1:end,:) = [UzP_top, UzUr_top, UzUz_top] - dt1/Re.*[UzP_top, UzUr_top, UzUz_top];
        b(end-n_magma+1:end) = (u1(end-n_magma+1:end)' + dt1/Re.*UzR_top);

    case 'BDF2'
        M(end-n_magma+1:end,:) = [UzP_top, UzUr_top, UzUz_top] - 1/Ft/Re.*[UzP_top, UzUr_top, UzUz_top];
        b(end-n_magma+1:end) = -Dt/Ft*u2(end-n_magma+1:end)' + Bt/Ft*u1(end-n_magma+1:end)' + 1/Ft/Re.*UzR_top;

     case 'Steady'
        M(end-n_magma+1:end,:) = -[UzP_top, UzUr_top, UzUz_top];
        b(end-n_magma+1:end) = UzR_top;
end

switch BC{4}
    case 'Dirichlet'
        M(end-n_magma+1:end,:) = 0;
        M(end-n_magma+1:end,end-n_magma+1:end) = eye(n_magma);
        b(end-n_magma+1:end) = BC{5}/U;
        
        M(end-(n_magma+1)*n_magma-n_magma:end-(n_magma+1)*n_magma,:) = 0;
        M(end-(n_magma+1)*n_magma-n_magma:end-(n_magma+1)*n_magma, ...
            end-(n_magma+1)*n_magma-n_magma:end-(n_magma+1)*n_magma) = eye(n_magma+1);
        b(end-(n_magma+1)*n_magma-n_magma:end-(n_magma+1)*n_magma) = 0;
end


% Center
if min(z_urr)<(max(z_urr)/(n_magma+1))
    M(length(P)+1:n_magma+1:length(P)+n_magma*(n_magma+1),:) = 0;
    M(length(P)+1:n_magma+1:length(P)+n_magma*(n_magma+1),length(P)+1:n_magma+1:length(P)+n_magma*(n_magma+1)) = eye(n_magma);
    b(length(P)+1:n_magma+1:length(P)+n_magma*(n_magma+1)) = 0;
else
    % Pressure at surface
    M(length(P) + 1:n_magma+1:(length(P) + n_magma*(n_magma+1)),:) = 0;
    
    dr = [z_pr(n_magma+1) - z_pr(1),...
        z_pr(2*n_magma+1:n_magma:end) - z_pr(n_magma+1:n_magma:end-2*n_magma+1), ...
        z_pr(end-n_magma+1) - z_pr(end-2*n_magma+1)];
    dz = [z_pz(n_magma+1) - z_pz(1),...
        z_pz(2*n_magma+1:n_magma:end) - z_pz(n_magma+1:n_magma:end-2*n_magma+1), ...
        z_pz(end-n_magma+1) - z_pz(end-2*n_magma+1)];
    thetas = atan2(dz,dr);
    
    UrP_left = zeros(n_magma,n_magma^2);
    UrP_left(:,1:n_magma:n_magma^2) = -1./rho_interpr(1:n_magma+1:end)'./(z_urr(1:n_magma+1:end)-z_p(1:n_magma:n_magma^2))'.*eye(n_magma);
    UrR_left = 1./rho_interpr(1:n_magma+1:end)'./(z_urr(1:n_magma+1:end)-z_p(1:n_magma:n_magma^2))'.*P0(1:n_magma:end)' - ...
        1./rho_interpr(1:n_magma+1:end)'.*(ddr_p(1:n_magma:n_magma^2,:)*P0' + 1./z_urr(1:n_magma+1:end)'.*P0(1:n_magma:end)');
    
    ind_left = find(ind_ur(1:(2*n_magma+1)^2));
    Ic = eye(n_magma*(n_magma+1));
    Iright2 = diag(ones(n_magma*(n_magma+1)-2,1),2);
    Iright1 = diag(ones(n_magma*(n_magma+1)-1,1),1);

    UrUr_left = 1./rho_interpr(1:n_magma+1:end)'.*1./z_urr(1:n_magma+1:end)'.*(diag(ddr_ur(1:n_magma+1:end,3:n_magma+1:end)).*(z_urr(3:n_magma+1:end)'.*etar(ind_left(3:n_magma+1:end))'.*(4/3*ddr_ur(3:n_magma+1:end,:) + ...
        (-2/3)./z_urr(3:n_magma+1:end)'.*Iright2(1:n_magma+1:end,:))) + ...
        diag(ddr_ur(1:n_magma+1:end,2:n_magma+1:end)).*(z_urr(2:n_magma+1:end)'.*etar(ind_left(2:n_magma+1:end))'.*(4/3*ddr_ur(2:n_magma+1:end,:) + ...
        (-2/3)./z_urr(2:n_magma+1:end)'.*Iright1(1:n_magma+1:end,:))) + ...
        diag(ddr_ur(1:n_magma+1:end,1:n_magma+1:end)).*(z_urr(1:n_magma+1:end)'.*etar(ind_left(1:n_magma+1:end))'.*tan(thetas-pi/2)'.*ddz_ur(1:n_magma+1:end,:))) + ...   
        -1./rho_interpr(1:n_magma+1:end)'./z_urr(1:n_magma+1:end)'.*etar(ind_left(1:n_magma+1:end))'.*(4/3./z_urr(1:n_magma+1:end)'.*Ic(1:n_magma+1:end,:) - 2/3*ddr_ur(1:n_magma+1:end,:)) + ...
        (1./rho_interpr(1:n_magma+1:end))'.*(etar(ind_left(1:n_magma+1:end))'.*d2dz_ur(1:n_magma+1:end,:) + ...
        detadz(ind_left(1:n_magma+1:end))'.*(ddz_ur(1:n_magma+1:end,:)));
    UrUz_left = 1./rho_interpr(1:n_magma+1:end)'.*1./z_urr(1:n_magma+1:end)'.*(diag(ddr_ur(1:n_magma+1:end,3:n_magma+1:end)).*(z_urr(3:n_magma+1:end)'.*etar(ind_left(3:n_magma+1:end))'.*(-2/3*duzdz_ur(3:n_magma+1:end,:))) + ...
        diag(ddr_ur(1:n_magma+1:end,2:n_magma+1:end)).*(z_urr(2:n_magma+1:end)'.*etar(ind_left(2:n_magma+1:end))'.*(-2/3*duzdz_ur(2:n_magma+1:end,:))) + ...
        diag(ddr_ur(1:n_magma+1:end,1:n_magma+1:end)).*(z_urr(1:n_magma+1:end)'.*etar(ind_left(1:n_magma+1:end))'.*tan(thetas-pi/2)'.*duzdr_ur(1:n_magma+1:end,:))) + ...   
        -1./rho_interpr(1:n_magma+1:end)'./z_urr(1:n_magma+1:end)'.*etar(ind_left(1:n_magma+1:end))'.*(-2/3*duzdz_ur(1:n_magma+1:end,:)) + ...
        (1./rho_interpr(1:n_magma+1:end))'.*(etar(ind_left(1:n_magma+1:end))'.*d2uzdrdz_ur(1:n_magma+1:end,:) + ...
        detadz(ind_left(1:n_magma+1:end))'.*(duzdr_ur(1:n_magma+1:end,:)));

    switch timescheme
        case 'BDF1'
            M(length(P) + 1:n_magma+1:(length(P) + n_magma*(n_magma+1)),:) = [UrP_left, UrUr_left, UrUz_left] - dt1/Re.*[UrP_left, UrUr_left, UrUz_left];
            b(length(P) + 1:n_magma+1:(length(P) + n_magma*(n_magma+1))) = (u1(length(P) + 1:n_magma+1:(length(P) + n_magma*(n_magma+1)))' + dt1/Re.*UrR_left);
    
        case 'BDF2'
            M(length(P) + 1:n_magma+1:(length(P) + n_magma*(n_magma+1)),:) = [UrP_left, UrUr_left, UrUz_left] - 1/Ft/Re.*[UrP_left, UrUr_left, UrUz_left];
            b(length(P) + 1:n_magma+1:(length(P) + n_magma*(n_magma+1))) = -Dt/Ft*u2(length(P) + 1:n_magma+1:(length(P) + n_magma*(n_magma+1)))' + Bt/Ft*u1(length(P) + 1:n_magma+1:(length(P) + n_magma*(n_magma+1)))' + 1/Ft/Re.*UrR_left;
    
         case 'Steady'
            M(length(P) + 1:n_magma+1:(length(P) + n_magma*(n_magma+1)),:) = -[UrP_left, UrUr_left, UrUz_left];
            b(length(P) + 1:n_magma+1:(length(P) + n_magma*(n_magma+1))) = UrR_left;
    end

    M(length(P)+1,:) = 0;
    M(length(P)+1,length(P)+1) = 1; 
    b(length(P)+1) = 0;

    M(length(P)+1+(n_magma+1)*(n_magma-1),:) = 0;
    M(length(P)+1+(n_magma+1)*(n_magma-1),length(P)+1+(n_magma+1)*(n_magma-1)) = 1; 
    b(length(P)+1+(n_magma+1)*(n_magma-1)) = 0;

end

%

% Right 
% Pressure at surface
M(length(P) + n_magma+1:n_magma+1:(length(P) + n_magma*(n_magma+1)),:) = 0;

dr = [z_pr(2*n_magma) - z_pr(n_magma),...
    z_pr(3*n_magma:n_magma:end) - z_pr(n_magma:n_magma:end-2*n_magma), ...
    z_pr(end) - z_pr(end-n_magma)];
dz = [z_pz(2*n_magma) - z_pz(n_magma),...
    z_pz(3*n_magma:n_magma:end) - z_pz(n_magma:n_magma:end-2*n_magma), ...
    z_pz(end) - z_pz(end-n_magma)];
thetas = atan2(dz,dr);

UrP_right = zeros(n_magma,n_magma^2);
UrP_right(:,n_magma:n_magma:n_magma^2) = -1./rho_interpr(n_magma+1:n_magma+1:end)'./(z_urr(n_magma+1:n_magma+1:end)-z_p(n_magma:n_magma:n_magma^2))'.*eye(n_magma);
UrR_right = 1./rho_interpr(n_magma+1:n_magma+1:end)'./(z_urr(n_magma+1:n_magma+1:end)-z_p(n_magma:n_magma:n_magma^2))'.*P0(n_magma:n_magma:end)' - ...
    1./rho_interpr(n_magma+1:n_magma+1:end)'.*(ddr_p(n_magma:n_magma:n_magma^2,:)*P0' + 1./z_urr(n_magma+1:n_magma+1:end)'.*P0(n_magma:n_magma:end)');

ind_right = find(ind_ur(1:(2*n_magma+1)^2));
Ic = eye(n_magma*(n_magma+1));
Ileft2 = diag(ones(n_magma*(n_magma+1)-2,1),-2);
Ileft1 = diag(ones(n_magma*(n_magma+1)-1,1),-1);

UrUr_right = 1./rho_interpr(n_magma+1:n_magma+1:end)'.*1./z_urr(n_magma+1:n_magma+1:end)'.*(diag(ddr_ur(n_magma+1:n_magma+1:end,n_magma-1:n_magma+1:end)).*(z_urr(n_magma-1:n_magma+1:end)'.*etar(ind_right(n_magma-1:n_magma+1:end))'.*(4/3*ddr_ur(n_magma-1:n_magma+1:end,:) + ...
    (-2/3)./z_urr(n_magma-1:n_magma+1:end)'.*Ileft2(n_magma+1:n_magma+1:end,:))) + ...
    diag(ddr_ur(n_magma+1:n_magma+1:end,n_magma:n_magma+1:end)).*(z_urr(n_magma:n_magma+1:end)'.*etar(ind_right(n_magma:n_magma+1:end))'.*(4/3*ddr_ur(n_magma:n_magma+1:end,:) + ...
    (-2/3)./z_urr(n_magma:n_magma+1:end)'.*Ileft1(n_magma+1:n_magma+1:end,:))) + ...
    diag(ddr_ur(n_magma+1:n_magma+1:end,n_magma+1:n_magma+1:end)).*(z_urr(n_magma+1:n_magma+1:end)'.*etar(ind_right(n_magma+1:n_magma+1:end))'.*tan(thetas-pi/2)'.*ddz_ur(n_magma+1:n_magma+1:end,:))) + ...   
    -1./rho_interpr(n_magma+1:n_magma+1:end)'./z_urr(n_magma+1:n_magma+1:end)'.*etar(ind_right(n_magma+1:n_magma+1:end))'.*(4/3./z_urr(n_magma+1:n_magma+1:end)'.*Ic(n_magma+1:n_magma+1:end,:) - 2/3*ddr_ur(n_magma+1:n_magma+1:end,:)) + ...
    (1./rho_interpr(n_magma+1:n_magma+1:end))'.*(etar(ind_right(n_magma+1:n_magma+1:end))'.*d2dz_ur(n_magma+1:n_magma+1:end,:) + ...
    detadz(ind_right(n_magma+1:n_magma+1:end))'.*(ddz_ur(n_magma+1:n_magma+1:end,:)));
UrUz_right = 1./rho_interpr(n_magma+1:n_magma+1:end)'.*1./z_urr(n_magma+1:n_magma+1:end)'.*(diag(ddr_ur(n_magma+1:n_magma+1:end,n_magma-1:n_magma+1:end)).*(z_urr(n_magma-1:n_magma+1:end)'.*etar(ind_right(n_magma-1:n_magma+1:end))'.*(-2/3*duzdz_ur(n_magma-1:n_magma+1:end,:))) + ...
    diag(ddr_ur(n_magma+1:n_magma+1:end,n_magma:n_magma+1:end)).*(z_urr(n_magma:n_magma+1:end)'.*etar(ind_right(n_magma:n_magma+1:end))'.*(-2/3*duzdz_ur(n_magma:n_magma+1:end,:))) + ...
    diag(ddr_ur(n_magma+1:n_magma+1:end,n_magma+1:n_magma+1:end)).*(z_urr(n_magma+1:n_magma+1:end)'.*etar(ind_right(n_magma+1:n_magma+1:end))'.*tan(thetas-pi/2)'.*duzdr_ur(n_magma+1:n_magma+1:end,:))) + ...   
    1./rho_interpr(n_magma+1:n_magma+1:end)'./z_urr(n_magma+1:n_magma+1:end)'.*etar(ind_right(n_magma+1:n_magma+1:end))'.*(2/3*duzdz_ur(n_magma+1:n_magma+1:end,:)) + ...
    (1./rho_interpr(n_magma+1:n_magma+1:end))'.*(etar(ind_right(n_magma+1:n_magma+1:end))'.*d2uzdrdz_ur(n_magma+1:n_magma+1:end,:) + ...
    detadz(ind_right(n_magma+1:n_magma+1:end))'.*(duzdr_ur(n_magma+1:n_magma+1:end,:)));

switch timescheme
    case 'BDF1'
        M(length(P) + n_magma+1:n_magma+1:(length(P) + n_magma*(n_magma+1)),:) = [UrP_right, UrUr_right, UrUz_right] - dt1/Re.*[UrP_right, UrUr_right, UrUz_right];
        b(length(P) + n_magma+1:n_magma+1:(length(P) + n_magma*(n_magma+1))) = (u1(length(P) + n_magma+1:n_magma+1:(length(P) + n_magma*(n_magma+1)))' + dt1/Re.*UrR_right);

    case 'BDF2'
        M(length(P) + n_magma+1:n_magma+1:(length(P) + n_magma*(n_magma+1)),:) = [UrP_right, UrUr_right, UrUz_right] - 1/Ft/Re.*[UrP_right, UrUr_right, UrUz_right];
        b(length(P) + n_magma+1:n_magma+1:(length(P) + n_magma*(n_magma+1))) = -Dt/Ft*u2(length(P) + n_magma+1:n_magma+1:(length(P) + n_magma*(n_magma+1)))' + Bt/Ft*u1(length(P) + n_magma+1:n_magma+1:(length(P) + n_magma*(n_magma+1)))' + 1/Ft/Re.*UrR_right;

     case 'Steady'
        M(length(P) + n_magma+1:n_magma+1:(length(P) + n_magma*(n_magma+1)),:) = -[UrP_right, UrUr_right, UrUz_right];
        b(length(P) + n_magma+1:n_magma+1:(length(P) + n_magma*(n_magma+1))) = UrR_right;
end

switch BC{2}
    case 'Confined'
        j = [n_magma+1, find(z_urr>=(R/L - 1e-5))];
        M(length(P) + j,:) = 0;
        M(length(P) + j,length(P) + j) = eye(length(j));
        b(length(P) + j) = 0;

        switch BC{3}
            case 'No slip'
                k = [n_magma floor(j(rem(j,n_magma+1)==0)/(n_magma+1))*n_magma+n_magma];

                M(length(P)+length(z_urr)+k,:) = 0;
                M(length(P)+length(z_urr)+k, ...
                  length(P)+length(z_urr)+k) = eye(length(k));
                b(length(P)+length(z_urr)+k) = 0;

                % Bottom right corner
                M(n_magma,:) = 0; 
                M(n_magma,n_magma) = -1;
                M(n_magma,n_magma-1) = 1/2;
                M(n_magma,2*n_magma) = 1/2;
                b(n_magma) = 0;

        end
end
Vi = V;
V = M\b;

P = V(1:length(P1))';
u = V(length(P1)+1:end)';

n = n+1;

if ~isreal(u) || ~isreal(P)
    %'STOP'
end

end

P = P/(L/Mu/U);
u = u*U;

%%%%%%%%%%%%% Spherical %%%%%%%%%%%%%%%%%
function [P,u] = Spherical_pde(P1,u1,P2,u2,P0,rho,drhodt,...
            phi,a,SurfTens,eta,beta,dbetadt,R,g,z_p,z_u,z_t,dt1,dt2,BC,timescheme)

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
z_t = z_t/L;
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

%z_t = reshape([z_u; z_p, 0],1,[]);
%z_t = z_t(1:end-1);

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

dudr = dudz*u' + 2./z_u'.*u'; 
dudr1 = dudz*u1' + 2./z_u'.*u1';
dudr2 = dudz*u2' + 2./z_u'.*u2';

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

% Pressure at surface
M(end,:) = 0;

%UP_end = 0*z_p;
%UP_end(end) = -1./rho_interp(end)./(z_u(end)-z_p(end));
%UR_end = 0*P0';
%UR_end = UR_end(end) + (1./rho_interp(end)).*(1./(z_u(end)-z_p(end)).*P0(end));
%UU_end = 0*z_u;
%UU_end(end-3:end) = 1./rho_interp(end)*4/3*[(z_u(end-2)./z_u(end)).^2.*D(end).*etar(end-2).*(-D(end-2)),...
%    (z_u(end-2)./z_u(end)).^2.*D(end)*etar(end-2).*(-E(end-2)-1./z_u(end-2)) + (z_u(end-1)./z_u(end)).^2.*(-B(end))*etar(end-1).*(-D(end-1)) + etar(end).*D(end),...
%    (z_u(end-2)./z_u(end)).^2.*D(end)*etar(end-2).*(C(end-2)) + (z_u(end-1)./z_u(end)).^2.*(-B(end))*etar(end-1).*(-E(end-1)-1./z_u(end-1)) + etar(end).*(-B(end)),...
%    (z_u(end-1)./z_u(end)).^2.*(-B(end))*etar(end-1).*C(end-1) + etar(end).*(F(end) - 1./z_u(end))];

UP_end = 0*z_p;
UP_end(end) = -1./rho_interp(end)./(z_u(end)-z_p(end));
UR_end = 0; %P0(end)./rho_interp(end)./(z_u(end)-z_p(end));
UU_end = 0*z_u;
eye0 = diag(ones(length(u)-1));
UU_end = 1./rho_interp(end).*(4/3.*(D(end)*etar(end-4)*dudz(end-2,:) + ...
    -B(end)*etar(end-2)*dudz(end-1,:) + ...
    -4*etar(end)/R^2*eye0(end,:)));

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
P = V(1:length(P))';
%if any(P<0)
%    P = P0
%end
u = V(length(P)+1:end)';

n = n+1;

end

P = P/(L/Mu/U);
u = u*U;

%%%%%%%%%%%%% 2D Planar %%%%%%%%%%%%%%%%%
function [P,u] = Planar2D_pde(P1,u1,P2,u2,P0,rho,drhodt,...
            phi,a,SurfTens,eta,beta,dbetadt,R,g,z_p,z_u,z_t,dt1,dt2,BC,timescheme)

% Non-dimensionalization
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
z_t = z_t/L;
dt1 = dt1/(L/U);
dt2 = dt2/(L/U);

% Time setpping scheme
Bt = (dt1 + dt2)/dt1/dt2;
Dt = dt1/dt2/(dt1+dt2);
Ft = (dt2+2*dt1)/dt1/(dt2+dt1);

% Initial conditions
P = P1; 
u = u1;
Vi = [P u]';
V = 0;

% Numerical tolerance
max_iter = 100;
tol = 1e-4;
n = 1;

% Split coordinates
n_magma = floor((sqrt(length(z_t)/2)-1)/2);

z_pr = z_p(1:n_magma^2);
z_pz = z_p(n_magma^2+1:end);

z_urr = z_u(1:n_magma*(n_magma+1));
z_urz = z_u(n_magma*(n_magma+1)+1:2*n_magma*(n_magma+1));
z_uzr = z_u(2*n_magma*(n_magma+1)+1:3*n_magma*(n_magma+1));
z_uzz = z_u(3*n_magma*(n_magma+1)+1:end);

% Calculate derivatives for variables
[ddr_p,ddz_p,d2dr_p,d2dz_p] = Derivs2D(z_p,n_magma,n_magma);
[ddr_t,ddz_t,d2dr_t,d2dz_t] = Derivs2D(z_t,2*n_magma+1,2*n_magma+1);
[ddr_ur,ddz_ur,d2dr_ur,d2dz_ur] = Derivs2D(z_u(1:2*n_magma*(n_magma+1)),n_magma+1,n_magma);
[ddr_uz,ddz_uz,d2dr_uz,d2dz_uz] = Derivs2D(z_u(2*n_magma*(n_magma+1)+1:end),n_magma,n_magma+1);

% Calculate derivative and interpolation of velocity on staggered nodes
% Radial
I = ones(size(z_urr));
I(n_magma+1:n_magma+1:end) = 0;

h1 = z_urr(2:end)-z_urr(1:end-1);
w1 = (z_pr - z_urr(logical(I)))./h1(logical(I(1:end-1)));
w2 = (z_urr(logical(circshift(I,1))) - z_pr)./h1(logical(I(1:end-1)));

ur_interp = zeros([n_magma^2,n_magma*(n_magma+1)]);
durdr_P = zeros([n_magma^2,n_magma*(n_magma+1)]);
for i = 1:n_magma
    durdr_P((i-1)*n_magma+1:i*n_magma,(i-1)*(n_magma+1)+1:i*(n_magma+1)-1) = -diag(1./h1((i-1)*(n_magma+1)+1:i*(n_magma+1)-1));
    durdr_P((i-1)*n_magma+1:i*n_magma,(i-1)*(n_magma+1)+2:i*(n_magma+1)) = durdr_P((i-1)*n_magma+1:i*n_magma,(i-1)*(n_magma+1)+2:i*(n_magma+1)) + ...
        diag(1./h1((i-1)*(n_magma+1)+1:i*(n_magma+1)-1));
    ur_interp((i-1)*n_magma+1:i*n_magma,(i-1)*(n_magma+1)+1:i*(n_magma+1)-1) = diag(w1((i-1)*n_magma+1:i*n_magma));
    ur_interp((i-1)*n_magma+1:i*n_magma,(i-1)*(n_magma+1)+2:i*(n_magma+1)) = ur_interp((i-1)*n_magma+1:i*n_magma,(i-1)*(n_magma+1)+2:i*(n_magma+1)) + ...
        diag(w2((i-1)*n_magma+1:i*n_magma));
end

% Vertical
z_uzz_trans = reshape(reshape(z_uzz,n_magma,n_magma+1)',1,[]);
h1 = z_uzz_trans(2:end)-z_uzz_trans(1:end-1);
h1 = reshape(reshape(h1(logical(I)),n_magma,n_magma)',1,[]);

w1 = reshape(reshape(reshape(reshape(z_pz,n_magma,n_magma)',1,[]) - z_uzz_trans(logical(I)),n_magma,n_magma)',1,[])./h1;
w2 = reshape(reshape(z_uzz_trans(logical(circshift(I,1))) - reshape(reshape(z_pz,n_magma,n_magma)',1,[]),n_magma,n_magma)',1,[])./h1;

uz_interp = zeros([n_magma^2,n_magma*(n_magma+1)]);
duzdz_P = zeros([n_magma^2,n_magma*(n_magma+1)]);
uz_interp(:,1:n_magma^2) = diag(w1);
uz_interp(:,n_magma+1:end) = uz_interp(:,n_magma+1:end) + diag(w2);
duzdz_P(:,1:n_magma^2) = diag(-1./h1);
duzdz_P(:,n_magma+1:end) = duzdz_P(:,n_magma+1:end) + diag(1./h1);

% Calculate derivatives on velocity nodes
% Radial

h1 = 0*z_urr;
h1(logical(circshift(I,1))) = [z_pr(2:end) - z_pr(1:end-1), 0];
h1(1:n_magma+1:end) = z_pr(1:n_magma:end);
h1(n_magma+1:n_magma+1:end) = z_urr(n_magma+1:n_magma+1:end) - z_pr(n_magma:n_magma:end);

dPdr_ur = zeros(n_magma*(n_magma+1),n_magma^2);
for i = 1:n_magma
    dPdr_ur((i-1)*(n_magma+1)+2:i*(n_magma+1)-1,(i-1)*(n_magma)+1:i*(n_magma)-1) = diag(-1./h1((i-1)*(n_magma+1)+2:i*(n_magma+1)-1));
    dPdr_ur((i-1)*(n_magma+1)+2:i*(n_magma+1)-1,(i-1)*(n_magma)+2:i*(n_magma)) = dPdr_ur((i-1)*(n_magma+1)+2:i*(n_magma+1)-1,(i-1)*(n_magma)+2:i*(n_magma)) + ...
        diag(1./h1((i-1)*(n_magma+1)+2:i*(n_magma+1)-1));
end
dPdr_ur(1:n_magma+1:end,:) = dPdr_ur(2:n_magma+1:end,:);
dPdr_ur(n_magma+1:n_magma+1:end,:) = dPdr_ur(n_magma:n_magma+1:end,:);

% Vertical

z_pz_trans = reshape(reshape(z_pz,n_magma,[])',1,[]);
h1 = reshape(reshape([z_pz(1) z_pz_trans(2:end) - z_pz_trans(1:end-1)],n_magma,[])',1,[]);
h1(1:n_magma) = z_pz(1:n_magma);
h1 = [h1 z_uzz(end-n_magma+1:end) - z_pz(end-n_magma+1:end)];

dPdz_uz = zeros(n_magma*(n_magma+1),n_magma^2);
dPdz_uz(n_magma+1:end-n_magma,1:end-n_magma) = diag(-1./h1(n_magma+1:end-n_magma)); 
dPdz_uz(n_magma+1:end-n_magma,n_magma+1:end) = dPdz_uz(n_magma+1:end-n_magma,n_magma+1:end) + diag(1./h1(n_magma+1:end-n_magma));

dPdz_uz(1:n_magma,:) = dPdz_uz(n_magma+1:2*n_magma,:);
dPdz_uz(end-n_magma+1:end,:) = dPdz_uz(end-2*n_magma+1:end-n_magma,:);

% Cross-derivatives in velocity
% Radial 

x1 = z_uzr(1:end-n_magma); y1 = z_uzz(1:end-n_magma);
x2 = circshift(z_uzr(1:end-n_magma),-1); y2 = z_uzz(2:end-n_magma+1);
x3 = circshift(z_uzr(n_magma+1:end),-1); y3 = circshift(z_uzz(n_magma+1:end),-1);
x4 = z_uzr(n_magma+1:end); y4 = z_uzz(n_magma+1:end);

x0 = z_urr(logical(circshift(I,1))); y0 = z_urz(logical(circshift(I,1)));

a_zr = (x1 - x2 + x3 - x4).*(-y1 + y4) - (-x1 + x4).*(y1 - y2 + y3 - y4);
b_zr = ((x1 - x2 + x3 - x4).*y1 - x1.*(y1 - y2 + y3 - y4) +...
    (-x1 + x2).*(-y1 + y4) - (-x1 + x4).*(-y1 + y2) + ...
    x0.*(y1 - y2 + y3 - y4) - y0.*(x1 - x2 + x3 - x4));
c_zr = (-x1 + x2).*y1 - x1.*(-y1 + y2) + x0.*(-y1 + y2) - y0.*(-x1 + x2);

m = (-b_zr + sqrt(b_zr.^2 - 4*a_zr.*c_zr))./(2*a_zr);
m(abs(a_zr)<1e-16) = -c_zr(abs(a_zr)<1e-16)./b_zr(abs(a_zr)<1e-16);
l = (x0 - x1 - (-x1 + x4).*m)./((-x1 + x2) + (x1 - x2 + x3 - x4).*m);

w1 = (z_urz(1:n_magma+1:end) - z_uzz(1:n_magma:end-n_magma-1))./(z_uzz(n_magma+1:n_magma:end) - z_uzz(1:n_magma:end-n_magma-1));
w2 = (z_uzz(n_magma+1:n_magma:end) - z_urz(1:n_magma+1:end))./(z_uzz(n_magma+1:n_magma:end) - z_uzz(1:n_magma:end-n_magma-1));
w3 = (z_urz(n_magma+1:n_magma+1:end) - z_uzz(n_magma:n_magma:end-n_magma))./(z_uzz(2*n_magma:n_magma:end) - z_uzz(n_magma:n_magma:end-n_magma));
w4 = (z_uzz(2*n_magma:n_magma:end) - z_urz(n_magma+1:n_magma+1:end))./(z_uzz(2*n_magma:n_magma:end) - z_uzz(n_magma:n_magma:end-n_magma));

uz_I = eye(length(z_uzr));
uz_interp_ur = zeros(n_magma*(n_magma+1));
uz_interp_ur(logical(circshift(I,1)),:) = ((1-l).*(1-m))'.*uz_I(1:end-n_magma,:) + ...
    ((1-l).*m)'.*uz_I(2:end-n_magma+1,:) + (l.*(1-m))'.*[uz_I(n_magma+2:end,:); uz_I(end,:)] + ...
    (l.*m)'.*uz_I(n_magma+1:end,:);

duzdr_ur = zeros(n_magma*(n_magma+1));
duzdr_ur(logical(circshift(I,1)),:) = ((1-l).*(1-m))'.*ddr_uz(1:end-n_magma,:) + ...
    ((1-l).*m)'.*ddr_uz(2:end-n_magma+1,:) + (l.*(1-m))'.*[ddr_uz(n_magma+2:end,:); ddr_uz(end,:)] + ...
    (l.*m)'.*ddr_uz(n_magma+1:end,:);

duzdz_ur = zeros(n_magma*(n_magma+1));
duzdz_ur(logical(circshift(I,1)),:) = ((1-l).*(1-m))'.*ddz_uz(1:end-n_magma,:) + ...
    ((1-l).*m)'.*ddz_uz(2:end-n_magma+1,:) + (l.*(1-m))'.*[ddz_uz(n_magma+2:end,:); ddz_uz(end,:)] + ...
    (l.*m)'.*ddz_uz(n_magma+1:end,:);

d2uzdz2_ur = zeros(n_magma*(n_magma+1));
d2uzdz2_ur(logical(circshift(I,1)),:) = ((1-l).*(1-m))'.*d2dz_uz(1:end-n_magma,:) + ...
    ((1-l).*m)'.*d2dz_uz(2:end-n_magma+1,:) + (l.*(1-m))'.*[d2dz_uz(n_magma+2:end,:); d2dz_uz(end,:)] + ...
    (l.*m)'.*d2dz_uz(n_magma+1:end,:);

uz_interp_ur(1:n_magma+1:end,:) = w1'.*uz_I(1:n_magma:end-n_magma-1,:) + w2'.*uz_I(n_magma+1:n_magma:end,:);
duzdr_ur(1:n_magma+1:end,:) = w1'.*ddr_uz(1:n_magma:end-n_magma-1,:) + w2'.*ddr_uz(n_magma+1:n_magma:end,:);
duzdr_ur(1:n_magma+1:end,:) = 0*duzdr_ur(1:n_magma+1:end,:);
duzdz_ur(1:n_magma+1:end,:) = w1'.*ddz_uz(1:n_magma:end-n_magma-1,:) + w2'.*ddz_uz(n_magma+1:n_magma:end,:);
d2uzdz2_ur(1:n_magma+1:end,:) = w1'.*d2dz_uz(1:n_magma:end-n_magma-1,:) + w2'.*d2dz_uz(n_magma+1:n_magma:end,:);

uz_interp_ur(n_magma+1:n_magma+1:end,:) = w3'.*uz_I(n_magma:n_magma:end-n_magma,:) + w4'.*uz_I(2*n_magma:n_magma:end,:);
duzdr_ur(n_magma+1:n_magma+1:end,:) = w3'.*ddr_uz(n_magma:n_magma:end-n_magma,:) + w4'.*ddr_uz(2*n_magma:n_magma:end,:);
duzdz_ur(n_magma+1:n_magma+1:end,:) = w3'.*ddz_uz(n_magma:n_magma:end-n_magma,:) + w4'.*ddz_uz(2*n_magma:n_magma:end,:);
d2uzdz2_ur(n_magma+1:n_magma+1:end,:) = w3'.*d2dz_uz(n_magma:n_magma:end-n_magma,:) + w4'.*d2dz_uz(2*n_magma:n_magma:end,:);
switch BC{2}
    case 'Confined'
        switch BC{3}
            case 'No slip'
                uz_interp_ur(n_magma+1:n_magma+1:end,:) = (w3'.*uz_I(n_magma:n_magma:end-n_magma,:) + w4'.*uz_I(2*n_magma:n_magma:end,:))/2;
                duzdr_ur(n_magma+1:n_magma+1:end,:) = (w3'.*ddr_uz(n_magma:n_magma:end-n_magma,:) + w4'.*ddr_uz(2*n_magma:n_magma:end,:))/2;
                duzdz_ur(n_magma+1:n_magma+1:end,:) = 0*duzdz_ur(n_magma+1:n_magma+1:end,:);
                d2uzdz2_ur(n_magma+1:n_magma+1:end,:) = 0*d2uzdz2_ur(n_magma+1:n_magma+1:end,:);
        end
end

d2uzdrdz_ur = ddz_ur*duzdr_ur;

% Vertical
x1 = z_urr(logical(I)); y1 = z_urz(logical(I));
x4 = x1(n_magma+1:end); y4 = y1(n_magma+1:end); 
x1 = x1(1:end-n_magma); y1 = y1(1:end-n_magma);

x2 = z_urr(logical(circshift(I,1))); y2 = z_urz(logical(circshift(I,1))); 
x3 = x2(n_magma+1:end); y3 = y2(n_magma+1:end); 
x2 = x2(1:end-n_magma); y2 = y2(1:end-n_magma);

x0 = z_uzr(n_magma+1:end-n_magma); y0 = z_uzz(n_magma+1:end-n_magma);

a_zz = (x1 - x2 + x3 - x4).*(-y1 + y4) - (-x1 + x4).*(y1 - y2 + y3 - y4);
b_zz = ((x1 - x2 + x3 - x4).*y1 - x1.*(y1 - y2 + y3 - y4) +...
    (-x1 + x2).*(-y1 + y4) - (-x1 + x4).*(-y1 + y2) + ...
    x0.*(y1 - y2 + y3 - y4) - y0.*(x1 - x2 + x3 - x4));
c_zz = (-x1 + x2).*y1 - x1.*(-y1 + y2) + x0.*(-y1 + y2) - y0.*(-x1 + x2);

m = (-b_zz + sqrt(b_zz.^2 - 4*a_zz.*c_zz))./(2*a_zz);
m(abs(a_zz)<1e-16) = -c_zz(abs(a_zz)<1e-16)./b_zz(abs(a_zz)<1e-16);
l = (x0 - x1 - (-x1 + x4).*m)./((-x1 + x2) + (x1 - x2 + x3 - x4).*m);
ur_I = eye(length(z_urr));

w1 = (z_uzr(1:n_magma) - z_urr(1:n_magma))./(z_urr(2:n_magma+1) - z_urr(1:n_magma));
w2 = (z_urr(2:n_magma+1) - z_uzr(1:n_magma))./(z_urr(2:n_magma+1) - z_urr(1:n_magma));
w3 = (z_uzr(end-n_magma+1:end) - z_urr(end-n_magma:end-1))./(z_urr(end-n_magma+1:end) - z_urr(end-n_magma:end-1));
w4 = (z_urr(end-n_magma+1:end) - z_uzr(end-n_magma+1:end))./(z_urr(end-n_magma+1:end) - z_urr(end-n_magma:end-1));

ur_interp_uz = zeros(n_magma*(n_magma+1));
durdr_uz = zeros(n_magma*(n_magma+1));
durdz_uz = zeros(n_magma*(n_magma+1));
for i = 1:n_magma-1
ur_interp_uz(i*n_magma+1:(i+1)*n_magma,:) = ((1-l((i-1)*n_magma+1:i*n_magma)).*(1-m((i-1)*n_magma+1:i*n_magma)))'.*ur_I((i-1)*(n_magma+1)+1:i*(n_magma+1)-1,:) + ...
    ((1-l((i-1)*n_magma+1:i*n_magma)).*m((i-1)*n_magma+1:i*n_magma))'.*ur_I((i-1)*(n_magma+1)+2:i*(n_magma+1),:) + ...
    (l((i-1)*n_magma+1:i*n_magma).*(1-m((i-1)*n_magma+1:i*n_magma)))'.*ur_I(i*(n_magma+1)+2:(i+1)*(n_magma+1),:) + ...
    (l((i-1)*n_magma+1:i*n_magma).*m((i-1)*n_magma+1:i*n_magma))'.*ur_I(i*(n_magma+1)+1:(i+1)*(n_magma+1)-1,:);

durdr_uz(i*n_magma+1:(i+1)*n_magma,:) = ((1-l((i-1)*n_magma+1:i*n_magma)).*(1-m((i-1)*n_magma+1:i*n_magma)))'.*ddr_ur((i-1)*(n_magma+1)+1:i*(n_magma+1)-1,:) + ...
    ((1-l((i-1)*n_magma+1:i*n_magma)).*m((i-1)*n_magma+1:i*n_magma))'.*ddr_ur((i-1)*(n_magma+1)+2:i*(n_magma+1),:) + ...
    (l((i-1)*n_magma+1:i*n_magma).*(1-m((i-1)*n_magma+1:i*n_magma)))'.*ddr_ur(i*(n_magma+1)+2:(i+1)*(n_magma+1),:) + ...
    (l((i-1)*n_magma+1:i*n_magma).*m((i-1)*n_magma+1:i*n_magma))'.*ddr_ur(i*(n_magma+1)+1:(i+1)*(n_magma+1)-1,:);

durdz_uz(i*n_magma+1:(i+1)*n_magma,:) = ((1-l((i-1)*n_magma+1:i*n_magma)).*(1-m((i-1)*n_magma+1:i*n_magma)))'.*ddz_ur((i-1)*(n_magma+1)+1:i*(n_magma+1)-1,:) + ...
    ((1-l((i-1)*n_magma+1:i*n_magma)).*m((i-1)*n_magma+1:i*n_magma))'.*ddz_ur((i-1)*(n_magma+1)+2:i*(n_magma+1),:) + ...
    (l((i-1)*n_magma+1:i*n_magma).*(1-m((i-1)*n_magma+1:i*n_magma)))'.*ddz_ur(i*(n_magma+1)+2:(i+1)*(n_magma+1),:) + ...
    (l((i-1)*n_magma+1:i*n_magma).*m((i-1)*n_magma+1:i*n_magma))'.*ddz_ur(i*(n_magma+1)+1:(i+1)*(n_magma+1)-1,:);
end

switch BC{1}
    case 'No stress'
        ur_interp_uz(1:n_magma,:) = w1'.*ur_I(1:n_magma,:) + w2'.*ur_I(2:n_magma+1,:);
    case 'No normal'
        ur_interp_uz(1:n_magma,:) = 0;
end
ur_interp_uz(end-n_magma+1:end,:) = w3'.*ur_I(end-n_magma:end-1,:) + w4'.*ur_I(end-n_magma+1:end,:);
switch BC{1}
    case 'No stress'
        durdr_uz(1:n_magma,:) = w1'.*ddr_ur(1:n_magma,:) + w2'.*ddr_ur(2:n_magma+1,:);
    case 'No normal'
        durdr_uz(1:n_magma,:) = 0;
end
durdr_uz(end-n_magma+1:end,:) = w3'.*ddr_ur(end-n_magma:end-1,:) + w4'.*ddr_ur(end-n_magma+1:end,:);

switch BC{1}
    case 'No stress'
        durdz_uz(1:n_magma,:) = w1'.*ddz_ur(1:n_magma,:) + w2'.*ddz_ur(2:n_magma+1,:);
        durdz_uz(end-n_magma+1:end,:) = w3'.*ddz_ur(end-n_magma:end-1,:) + w4'.*ddz_ur(end-n_magma+1:end,:);
    case 'No normal'
        durdz_uz(1:n_magma,:) = (w1'.*ddz_ur(1:n_magma,:) + w2'.*ddz_ur(2:n_magma+1,:))/2;
        durdz_uz(end-n_magma+1:end,:) = (w3'.*ddz_ur(end-n_magma:end-1,:) + w4'.*ddz_ur(end-n_magma+1:end,:))/2;
end

d2urdrdz_uz = ddr_uz*durdz_uz;

% Calculate derivatives and interpolations of constants
dP0dz = (ddz_p*P0')';
drhodr = (ddr_t*rho')';
drhodz = (ddz_t*rho')';
dbetadr = (ddz_p*beta')';
dbetadz = (ddz_p*beta')';

a_interp = scatteredInterpolant(z_pr',z_pz',a','linear','nearest');
a_interp = a_interp(z_t(1:(2*n_magma+1)^2)', z_t((2*n_magma+1)^2+1:end)')';
phi_interp = scatteredInterpolant(z_pr',z_pz',phi','linear','nearest');
phi_interp = phi_interp(z_t(1:(2*n_magma+1)^2)', z_t((2*n_magma+1)^2+1:end)')';

z_t_r = reshape(z_t(1:(2*n_magma+1)^2),2*n_magma+1,2*n_magma+1);
z_t_z = reshape(z_t((2*n_magma+1)^2+1:end),2*n_magma+1,2*n_magma+1);
ind_p = 0*z_t_r;
ind_p(2:2:end,2:2:end) = 1; ind_p = logical(repmat(reshape(ind_p,1,[]),1,2));
ind_ur = 0*z_t_r;
ind_ur(1:2:end,2:2:end) = 1; ind_ur = logical(repmat(reshape(ind_ur,1,[]),1,2));
ind_uz = 0*z_t_z;
ind_uz(2:2:end,1:2:end) = 1; ind_uz = logical(repmat(reshape(ind_uz,1,[]),1,2));

eta_interp = eta;
rho_interpr = rho(ind_ur(1:(2*n_magma+1)^2));
rho_interpz = rho(ind_uz(1:(2*n_magma+1)^2));

while (n<max_iter && norm(V-Vi)>tol) || n==1

PP = diag(-1./rho(ind_p(1:(2*n_magma+1)^2)).*(drhodt + ...
     1/2.*drhodr(ind_p(1:(2*n_magma+1)^2)).*(ur_interp*u(1:n_magma*(n_magma+1))')' + ...
     1/2.*drhodz(ind_p(1:(2*n_magma+1)^2)).*(uz_interp*u(n_magma*(n_magma+1)+1:end)')')) + ...
     diag(-1./beta.*(dbetadt + 1/2.*dbetadr.*(ur_interp*u(1:n_magma*(n_magma+1))')' + ...
     1/2.*dbetadz.*(uz_interp*u(n_magma*(n_magma+1)+1:end)')')) + ...
     (-1/2*(ur_interp*u(1:n_magma*(n_magma+1))').*ddr_p) + ...
     (-1/2*(uz_interp*u(n_magma*(n_magma+1)+1:end)').*(ddz_p));
PR = (-1./rho(ind_p(1:(2*n_magma+1)^2)).*(drhodt).*(1./beta-P0))' +...
    ((-1./beta).*dbetadt.*(-P0))';
PUr = (-1./rho(ind_p(1:(2*n_magma+1)^2)).*drhodr(ind_p(1:(2*n_magma+1)^2)).*(1./beta + 1/2*P-P0))'.*ur_interp + ...
    (-1./beta.*dbetadr.*(1/2*P-P0))'.*ur_interp + ...
    (-1/2*ddr_p*P').*ur_interp + ...
    (-(1./beta + 1/2*P-P0))'.*(durdr_P);
PUz = (-1./rho(ind_p(1:(2*n_magma+1)^2)).*drhodz(ind_p(1:(2*n_magma+1)^2)).*(1./beta + 1/2*P-P0))'.*uz_interp + ...
    (-1./beta.*dbetadz.*(1/2*P-P0))'.*uz_interp + ...
    (-1/2*(ddz_p*P') + dP0dz').*uz_interp + ...
    (-(1./beta + 1/2*P-P0))'.*duzdz_P;

% Calculate effective viscosity

% rinv = (1./z_urr);
% rinv(~isfinite(1./z_urr)) = 0;
% 
% gamma_ur = ddr_ur*u(1:n_magma*(n_magma+1))' + 1/2*ddz_ur*u(1:n_magma*(n_magma+1))' +...
%            (rinv.*u(1:n_magma*(n_magma+1)))' + 1/2*duzdr_ur*u(n_magma*(n_magma+1)+1:end)' + ...
%            duzdz_ur*u(n_magma*(n_magma+1)+1:end)';
% gamma_ur1 = ddr_ur*u1(1:n_magma*(n_magma+1))' + 1/2*ddz_ur*u1(1:n_magma*(n_magma+1))' +...
%            (rinv.*u1(1:n_magma*(n_magma+1)))' + 1/2*duzdr_ur*u1(n_magma*(n_magma+1)+1:end)' + ...
%            duzdz_ur*u1(n_magma*(n_magma+1)+1:end)';
% gamma_ur2 = ddr_ur*u2(1:n_magma*(n_magma+1))' + 1/2*ddz_ur*u2(1:n_magma*(n_magma+1))' +...
%            (rinv.*u2(1:n_magma*(n_magma+1)))' + 1/2*duzdr_ur*u2(n_magma*(n_magma+1)+1:end)' + ...
%            duzdz_ur*u2(n_magma*(n_magma+1)+1:end)';
% 
% gamma_uz = durdr_uz*u(1:n_magma*(n_magma+1))' + 1/2*durdz_uz*u(1:n_magma*(n_magma+1))' +...
%            (1./z_uzr)'.*(ur_interp_uz*u(1:n_magma*(n_magma+1))') + 1/2*ddr_uz*u(n_magma*(n_magma+1)+1:end)' + ...
%            ddz_uz*u(n_magma*(n_magma+1)+1:end)';
% gamma_uz1 = durdr_uz*u1(1:n_magma*(n_magma+1))' + 1/2*durdz_uz*u1(1:n_magma*(n_magma+1))' +...
%            (1./z_uzr)'.*(ur_interp_uz*u1(1:n_magma*(n_magma+1))') + 1/2*ddr_uz*u1(n_magma*(n_magma+1)+1:end)' + ...
%            ddz_uz*u1(n_magma*(n_magma+1)+1:end)';
% gamma_uz2 = durdr_uz*u2(1:n_magma*(n_magma+1))' + 1/2*durdz_uz*u2(1:n_magma*(n_magma+1))' +...
%            (1./z_uzr)'.*(ur_interp_uz*u2(1:n_magma*(n_magma+1))') + 1/2*ddr_uz*u2(n_magma*(n_magma+1)+1:end)' + ...
%            ddz_uz*u2(n_magma*(n_magma+1)+1:end)';

gamma_ur = ddz_ur*u(1:n_magma*(n_magma+1))' + uz_interp_ur*(ddr_uz*u(n_magma*(n_magma+1)+1:end)');
gamma_uz = ur_interp_uz*(ddz_ur*u(1:n_magma*(n_magma+1))') + ddr_uz*u(n_magma*(n_magma+1)+1:end)';

gamma_ur1 = ddz_ur*u1(1:n_magma*(n_magma+1))' + uz_interp_ur*(ddr_uz*u1(n_magma*(n_magma+1)+1:end)');
gamma_uz1 = ur_interp_uz*(ddz_ur*u1(1:n_magma*(n_magma+1))') + ddr_uz*u1(n_magma*(n_magma+1)+1:end)';

gamma_ur2 = ddz_ur*u2(1:n_magma*(n_magma+1))' + uz_interp_ur*(ddr_uz*u2(n_magma*(n_magma+1)+1:end)');
gamma_uz2 = ur_interp_uz*(ddz_ur*u2(1:n_magma*(n_magma+1))') + ddr_uz*u2(n_magma*(n_magma+1)+1:end)';

switch timescheme
    case 'BDF1'
        d2urdrdt = 0*gamma_ur;
        d2uzdrdt = 0*gamma_uz;
    case 'BDF2'
        d2urdrdt = Ft*gamma_ur -Bt*gamma_ur1 + Dt*gamma_ur2;
        d2uzdrdt = Ft*gamma_uz -Bt*gamma_uz1 + Dt*gamma_uz2;
    case 'Steady'
        d2urdrdt = 0;
        d2uzdrdt = 0;
end

gamma_interp = scatteredInterpolant([z_urr,z_uzr]',[z_urz,z_uzz]',[gamma_ur;gamma_uz],'linear','nearest');
gamma_interp = gamma_interp(z_t(1:(2*n_magma+1)^2),z_t((2*n_magma+1)^2+1:end));
d2udrdt_interp = scatteredInterpolant([z_urr,z_uzr]',[z_urz,z_uzz]',[d2urdrdt;d2uzdrdt],'linear','nearest');
d2udrdt_interp = d2udrdt_interp(z_t(1:(2*n_magma+1)^2),z_t((2*n_magma+1)^2+1:end));

Cc = max(sqrt((d2udrdt_interp'./gamma_interp').^2 + gamma_interp'.^2)./(L/U),1e-10)'.*a_interp'.*eta_interp*Mu/SurfTens;
eta0 = (1-phi_interp').^(-1);
etainf = (1-phi_interp').^(5/3);
etar = eta_interp.*(etainf + (eta0-etainf)./(1+(6/5*Cc).^2));
etar(etar>1e12/Mu) = 1e12/Mu;

detadr = ddr_t*etar';
detadz = ddz_t*etar';

UrP = -(1./rho_interpr').*dPdr_ur;
UrR = 0*z_urr';
UrUr = 4/3*(1./rho_interpr)'.*(etar(ind_ur(1:(2*n_magma+1)^2))'.*d2dr_ur + ...
    detadr(ind_ur(1:(2*n_magma+1)^2))'.*ddr_ur) + ...
    (1./rho_interpr)'.*etar(ind_ur(1:(2*n_magma+1)^2))'.*(d2dz_ur) + ...
    (1./rho_interpr)'.*detadz(ind_ur(1:(2*n_magma+1)^2))'.*(ddz_ur);
UrUz = (1./rho_interpr)'.*(etar(ind_ur(1:(2*n_magma+1)^2))'.*(1/3*d2uzdrdz_ur) + ...
    detadr(ind_ur(1:(2*n_magma+1)^2))'.*(-2/3*duzdz_ur) + ...
    detadz(ind_ur(1:(2*n_magma+1)^2))'.*(duzdr_ur));

UzP = -(1./rho_interpz').*dPdz_uz;
UzR = (1./rho_interpz').*0; %(2*rho_interpz'*g - dPdz_uz*P0');
UzUr = (1./rho_interpz)'.*(etar(ind_uz(1:(2*n_magma+1)^2))'.*(1/3*d2urdrdz_uz) + ...
        detadr(ind_uz(1:(2*n_magma+1)^2))'.*durdz_uz + ...
        detadz(ind_uz(1:(2*n_magma+1)^2))'.*(-2/3*durdr_uz));
UzUz = (1./rho_interpz)'.*(etar(ind_uz(1:(2*n_magma+1)^2))'.*(4/3*d2dz_uz + ...
        d2dr_uz) + ...
        detadr(ind_uz(1:(2*n_magma+1)^2))'.*ddr_uz + ...
        detadz(ind_uz(1:(2*n_magma+1)^2))'.*(4/3*ddz_uz));

switch timescheme
    case 'BDF1'
    M = eye(size([PP, PUr, PUz; UrP, UrUr, UrUz; UzP, UzUr, UzUz])) - dt1/Re.*[PP, PUr, PUz; UrP, UrUr, UrUz; UzP, UzUr, UzUz];
    b = ([P1 u1]' + dt1/Re.*[PR; UrR; UzR]);

    case 'BDF2'
        M = eye(size([PP, PUr, PUz; UrP, UrUr, UrUz; UzP, UzUr, UzUz])) - 1/Ft/Re.*[PP, PUr, PUz; UrP, UrUr, UrUz; UzP, UzUr, UzUz];
        b = (-Dt/Ft*[P2 u2]' + Bt/Ft*[P1 u1]' + 1/Ft/Re.*[PR; UrR; UzR]);

    case 'Steady'
        M = -[PP, PUr, PUz; UrP, UrUr, UrUz; UzP, UzUr, UzUz];
        b = [PR; UrR; UzR];
end

% BCs
% Velocity at bottom
switch BC{1}
    case 'No normal'
        M(length(P)+1:length(P)+n_magma+1,:) = 0;
        M(length(P)+1:length(P)+n_magma+1,length(P)+1:length(P)+n_magma+1) = eye(n_magma+1);
        b(length(P)+1:length(P)+n_magma+1) = 0;

        M(length(P)+length(z_urr)+1:length(P)+length(z_urr)+n_magma,:) = 0;
        M(length(P)+length(z_urr)+1:length(P)+length(z_urr)+n_magma,...
            length(P)+length(z_urr)+1:length(P)+length(z_urr)+n_magma) = eye(n_magma);
        b(length(P)+length(z_urr)+1:length(P)+length(z_urr)+n_magma) = 0;

    case 'No stress'
        %M(length(P)+1:length(P)+n_magma+1,:) = 0;
        %M(length(P)+1:length(P)+n_magma+1,length(P)+1:length(P)+n_magma+1) = eye(n_magma+1);
        %b(length(P)+1:length(P)+n_magma+1) = 0;

        M(length(P)+length(z_urr)+1:length(P)+length(z_urr)+n_magma,:) = 0;
        M(length(P)+length(z_urr)+1:length(P)+length(z_urr)+n_magma,...
            length(P)+length(z_urr)+1:length(P)+length(z_urr)+n_magma) = -eye(n_magma);
        M(length(P)+length(z_urr)+1:length(P)+length(z_urr)+n_magma,...
            length(P)+length(z_urr)+n_magma+1:length(P)+length(z_urr)+2*n_magma) = eye(n_magma);
        b(length(P)+length(z_urr)+1:length(P)+length(z_urr)+n_magma) = 0;
end

% Pressure at surface
M(end-n_magma+1:end,:) = 0;

dr = [z_pr(end-n_magma+2) - z_pr(end-n_magma+1),...
    z_pr(end-n_magma+3:end) - z_pr(end-n_magma+1:end-2), ...
    z_pr(end) - z_pr(end-1)];
dz = [z_pz(end-n_magma+2) - z_pz(end-n_magma+1),...
    z_pz(end-n_magma+3:end) - z_pz(end-n_magma+1:end-2), ...
    z_pz(end) - z_pz(end-1)];
thetas = atan(dz./dr);

UzP_top = zeros(n_magma,n_magma.^2);
UzP_top(:,end-n_magma+1:end) = -1./rho_interpz(end-n_magma+1:end)'./(z_uzz(end-n_magma+1:end)-z_p(end-n_magma+1:end))'.*eye(n_magma);
UzR_top = 1./rho_interpz(end-n_magma+1:end)'./(z_uzz(end-n_magma+1:end)-z_p(end-n_magma+1:end))'.*P0(end-n_magma+1:end)'; % - ...
    %2./rho_interpz(end-n_magma+1:end)'.*(ddz_p(end-n_magma+1:end,:)*P0') + 2*g;

ind_top = find(ind_uz(1:(2*n_magma+1)^2));

UzUr_top = 1./rho_interpz(end-n_magma+1:end)'*2/3.*(diag(ddz_uz(end-n_magma+1:end,end-3*n_magma+1:end-n_magma)).*(etar(ind_top(end-3*n_magma+1:end-2*n_magma))'.*(-durdr_uz(end-3*n_magma+1:end-2*n_magma,:))) + ...
    diag(ddz_uz(end-n_magma+1:end,end-2*n_magma+1:end-n_magma)).*(etar(ind_top(end-2*n_magma+1:end-n_magma))'.*(-durdr_uz(end-2*n_magma+1:end-n_magma,:)))) + ...
    1./rho_interpz(end-n_magma+1:end)'.*diag(ddz_uz(end-n_magma+1:end,end-n_magma+1:end)).*(etar(ind_top(end-n_magma+1:end))'.*tan(thetas)'.*durdz_uz(end-n_magma+1:end,:)) + ...
    1./rho_interpz(end-n_magma+1:end)'.*etar(ind_top(end-n_magma+1:end))'.*(d2urdrdz_uz(end-n_magma+1:end,:)) + ...
    1./rho_interpz(end-n_magma+1:end)'.*detadr(ind_top(end-n_magma+1:end))'.*durdz_uz(end-n_magma+1:end,:);

UzUz_top = 1./rho_interpz(end-n_magma+1:end)'*4/3.*(diag(ddz_uz(end-n_magma+1:end,end-3*n_magma+1:end-2*n_magma)).*(etar(ind_top(end-3*n_magma+1:end-2*n_magma))'.*(ddz_uz(end-3*n_magma+1:end-2*n_magma,:))) + ...
    diag(ddz_uz(end-n_magma+1:end,end-2*n_magma+1:end)).*(etar(ind_top(end-2*n_magma+1:end-n_magma))'.*(ddz_uz(end-2*n_magma+1:end-n_magma,:)))) + ...
    1./rho_interpz(end-n_magma+1:end)'.*diag(ddz_uz(end-n_magma+1:end,end-n_magma+1:end)).*(etar(ind_top(end-n_magma+1:end))'.*tan(thetas)'.*ddr_uz(end-n_magma+1:end,:)) + ...
    1./rho_interpz(end-n_magma+1:end)'.*etar(ind_top(end-n_magma+1:end))'.*(1./z_uzr(end-n_magma+1:end)'.*ddr_uz(end-n_magma+1:end,:) + d2dr_uz(end-n_magma+1:end,:)) + ...
    1./rho_interpz(end-n_magma+1:end)'.*detadr(ind_top(end-n_magma+1:end))'.*ddr_uz(end-n_magma+1:end,:);


switch timescheme
    case 'BDF1'
        M(end-n_magma+1:end,:) = [UzP_top, UzUr_top, UzUz_top] - dt1/Re.*[UzP_top, UzUr_top, UzUz_top];
        b(end-n_magma+1:end) = (u1(end-n_magma+1:end)' + dt1/Re.*UzR_top);

    case 'BDF2'
        M(end-n_magma+1:end,:) = [UzP_top, UzUr_top, UzUz_top] - 1/Ft/Re.*[UzP_top, UzUr_top, UzUz_top];
        b(end-n_magma+1:end) = -Dt/Ft*u2(end-n_magma+1:end)' + Bt/Ft*u1(end-n_magma+1:end)' + 1/Ft/Re.*UzR_top;

     case 'Steady'
        M(end-n_magma+1:end,:) = -[UzP_top, UzUr_top, UzUz_top];
        b(end-n_magma+1:end) = UzR_top;
end

%M(length(P)-n_magma+1,:) = 0;
%M(length(P)-n_magma+1,length(P)-n_magma+1) = 1;
%b(length(P)-n_magma+1) = P0(length(P)-n_magma+1);

% M(length(P)-n_magma+1:length(P),:) = 0;
% M(length(P)-n_magma+1:length(P),length(P)-n_magma+1:length(P)) = eye(n_magma);
% b(length(P)-n_magma+1:length(P)) = P0(length(P)-n_magma+1:end);

% Center
M(length(P)+1:n_magma+1:length(P)+n_magma*(n_magma+1),:) = 0;
M(length(P)+1:n_magma+1:length(P)+n_magma*(n_magma+1),length(P)+1:n_magma+1:length(P)+n_magma*(n_magma+1)) = eye(n_magma);
b(length(P)+1:n_magma+1:length(P)+n_magma*(n_magma+1)) = 0;

%M(length(P)+length(z_urr)+1:n_magma:end,:) = 0;
%M(length(P)+length(z_urr)+1:n_magma:end, ...
%  length(P)+length(z_urr)+1:n_magma:end) = eye(n_magma+1);
%b(length(P)+length(z_urr)+1:n_magma:end) = 0;


% Right 
% Pressure at surface
M(length(P) + n_magma+1:n_magma+1:(length(P) + n_magma*(n_magma+1)),:) = 0;

dr = [z_pr(2*n_magma) - z_pr(n_magma),...
    z_pr(3*n_magma:n_magma:end) - z_pr(n_magma:n_magma:end-2*n_magma), ...
    z_pr(end) - z_pr(end-n_magma)];
dz = [z_pz(2*n_magma) - z_pz(n_magma),...
    z_pz(3*n_magma:n_magma:end) - z_pz(n_magma:n_magma:end-2*n_magma), ...
    z_pz(end) - z_pz(end-n_magma)];
thetas = atan(dz./dr);

UrP_right = zeros(n_magma,n_magma^2);
UrP_right(:,n_magma:n_magma:n_magma^2) = -1./rho_interpr(n_magma+1:n_magma+1:end)'./(z_urr(n_magma+1:n_magma+1:end)-z_p(n_magma:n_magma:n_magma^2))'.*eye(n_magma);
UrR_right = 1./rho_interpr(n_magma+1:n_magma+1:end)'./(z_urr(n_magma+1:n_magma+1:end)-z_p(n_magma:n_magma:n_magma^2))'.*P0(n_magma:n_magma:end)' - ...
    1./rho_interpr(n_magma+1:n_magma+1:end)'.*(ddr_p(n_magma:n_magma:n_magma^2,:)*P0' + 1./z_urr(n_magma+1:n_magma+1:end)'.*P0(n_magma:n_magma:end)');

ind_right = find(ind_ur(1:(2*n_magma+1)^2));
Ic = eye(n_magma*(n_magma+1));
Ileft2 = diag(ones(n_magma*(n_magma+1)-2,1),-2);
Ileft1 = diag(ones(n_magma*(n_magma+1)-1,1),-1);

% UzUr_top = 1./rho_interpz(end-n_magma+1:end)'*2/3.*(diag(ddz_uz(end-n_magma+1:end,end-3*n_magma+1:end-n_magma)).*(etar(ind_top(end-3*n_magma+1:end-2*n_magma))'.*(-durdr_uz(end-3*n_magma+1:end-2*n_magma,:) - ...
%     1./z_uzr(end-3*n_magma+1:end-2*n_magma)'.*ur_interp_uz(end-3*n_magma+1:end-2*n_magma,:))) + ...
%     diag(ddz_uz(end-n_magma+1:end,end-2*n_magma+1:end-n_magma)).*(etar(ind_top(end-2*n_magma+1:end-n_magma))'.*(-durdr_uz(end-2*n_magma+1:end-n_magma,:)) - ...
%     1./z_uzr(end-2*n_magma+1:end-n_magma)'.*ur_interp_uz(end-2*n_magma+1:end-n_magma,:))) + ...
%     1./rho_interpz(end-n_magma+1:end)'.*diag(ddz_uz(end-n_magma+1:end,end-n_magma+1:end)).*(etar(ind_top(end-n_magma+1:end))'.*tan(thetas)'.*durdz_uz(end-n_magma+1:end,:)) + ...
%     1./rho_interpz(end-n_magma+1:end)'.*etar(ind_top(end-n_magma+1:end))'.*(1./z_uzr(end-n_magma+1:end)'.*durdz_uz(end-n_magma+1:end,:) + d2urdrdz_uz(end-n_magma+1:end,:)) + ...
%     1./rho_interpz(end-n_magma+1:end)'.*detadr(ind_top(end-n_magma+1:end))'.*durdz_uz(end-n_magma+1:end,:);
% 
% UzUz_top = 1./rho_interpz(end-n_magma+1:end)'*4/3.*(diag(ddz_uz(end-n_magma+1:end,end-3*n_magma+1:end-2*n_magma)).*(etar(ind_top(end-3*n_magma+1:end-2*n_magma))'.*(ddz_uz(end-3*n_magma+1:end-2*n_magma,:))) + ...
%     diag(ddz_uz(end-n_magma+1:end,end-2*n_magma+1:end)).*(etar(ind_top(end-2*n_magma+1:end-n_magma))'.*(ddz_uz(end-2*n_magma+1:end-n_magma,:)))) + ...
%     1./rho_interpz(end-n_magma+1:end)'.*diag(ddz_uz(end-n_magma+1:end,end-n_magma+1:end)).*(etar(ind_top(end-n_magma+1:end))'.*tan(thetas)'.*ddr_uz(end-n_magma+1:end,:)) + ...
%     1./rho_interpz(end-n_magma+1:end)'.*etar(ind_top(end-n_magma+1:end))'.*(1./z_uzr(end-n_magma+1:end)'.*ddr_uz(end-n_magma+1:end,:) + d2dr_uz(end-n_magma+1:end,:)) + ...
%     1./rho_interpz(end-n_magma+1:end)'.*detadr(ind_top(end-n_magma+1:end))'.*ddr_uz(end-n_magma+1:end,:);


UrUr_right = 1./rho_interpr(n_magma+1:n_magma+1:end)'.*1./z_urr(n_magma+1:n_magma+1:end)'.*(diag(ddr_ur(n_magma+1:n_magma+1:end,n_magma-1:n_magma+1:end)).*(z_urr(n_magma-1:n_magma+1:end)'.*etar(ind_right(n_magma-1:n_magma+1:end))'.*(4/3*ddr_ur(n_magma-1:n_magma+1:end,:) + ...
    (-2/3)./z_urr(n_magma-1:n_magma+1:end)'.*Ileft2(n_magma+1:n_magma+1:end,:))) + ...
    diag(ddr_ur(n_magma+1:n_magma+1:end,n_magma:n_magma+1:end)).*(z_urr(n_magma:n_magma+1:end)'.*etar(ind_right(n_magma:n_magma+1:end))'.*(4/3*ddr_ur(n_magma:n_magma+1:end,:) + ...
    (-2/3)./z_urr(n_magma:n_magma+1:end)'.*Ileft1(n_magma+1:n_magma+1:end,:))) + ...
    diag(ddr_ur(n_magma+1:n_magma+1:end,n_magma+1:n_magma+1:end)).*(z_urr(n_magma+1:n_magma+1:end)'.*etar(ind_right(n_magma+1:n_magma+1:end))'.*1./tan(thetas)'.*ddz_ur(n_magma+1:n_magma+1:end,:))) + ...   
    -1./rho_interpr(n_magma+1:n_magma+1:end)'./z_urr(n_magma+1:n_magma+1:end)'.*etar(ind_right(n_magma+1:n_magma+1:end))'.*(4/3./z_urr(n_magma+1:n_magma+1:end)'.*Ic(n_magma+1:n_magma+1:end,:) - 2/3*ddr_ur(n_magma+1:n_magma+1:end,:)) + ...
    (1./rho_interpr(n_magma+1:n_magma+1:end))'.*(etar(ind_right(n_magma+1:n_magma+1:end))'.*d2dz_ur(n_magma+1:n_magma+1:end,:) + ...
    detadz(ind_right(n_magma+1:n_magma+1:end))'.*(ddz_ur(n_magma+1:n_magma+1:end,:)));
UrUz_right = 1./rho_interpr(n_magma+1:n_magma+1:end)'.*1./z_urr(n_magma+1:n_magma+1:end)'.*(diag(ddr_ur(n_magma+1:n_magma+1:end,n_magma-1:n_magma+1:end)).*(z_urr(n_magma-1:n_magma+1:end)'.*etar(ind_right(n_magma-1:n_magma+1:end))'.*(-2/3*duzdz_ur(n_magma-1:n_magma+1:end,:))) + ...
    diag(ddr_ur(n_magma+1:n_magma+1:end,n_magma:n_magma+1:end)).*(z_urr(n_magma:n_magma+1:end)'.*etar(ind_right(n_magma:n_magma+1:end))'.*(-2/3*duzdz_ur(n_magma:n_magma+1:end,:))) + ...
    diag(ddr_ur(n_magma+1:n_magma+1:end,n_magma+1:n_magma+1:end)).*(z_urr(n_magma+1:n_magma+1:end)'.*etar(ind_right(n_magma+1:n_magma+1:end))'.*1./tan(thetas)'.*duzdr_ur(n_magma+1:n_magma+1:end,:))) + ...   
    -1./rho_interpr(n_magma+1:n_magma+1:end)'./z_urr(n_magma+1:n_magma+1:end)'.*etar(ind_right(n_magma+1:n_magma+1:end))'.*(-2/3*duzdz_ur(n_magma+1:n_magma+1:end,:)) + ...
    (1./rho_interpr(n_magma+1:n_magma+1:end))'.*(etar(ind_right(n_magma+1:n_magma+1:end))'.*d2uzdrdz_ur(n_magma+1:n_magma+1:end,:) + ...
    detadz(ind_right(n_magma+1:n_magma+1:end))'.*(duzdr_ur(n_magma+1:n_magma+1:end,:)));

switch timescheme
    case 'BDF1'
        M(length(P) + n_magma+1:n_magma+1:(length(P) + n_magma*(n_magma+1)),:) = [UrP_right, UrUr_right, UrUz_right] - dt1/Re.*[UrP_right, UrUr_right, UrUz_right];
        b(length(P) + n_magma+1:n_magma+1:(length(P) + n_magma*(n_magma+1))) = (u1(length(P) + n_magma+1:n_magma+1:(length(P) + n_magma*(n_magma+1)))' + dt1/Re.*UrR_right);

    case 'BDF2'
        M(length(P) + n_magma+1:n_magma+1:(length(P) + n_magma*(n_magma+1)),:) = [UrP_right, UrUr_right, UrUz_right] - 1/Ft/Re.*[UrP_right, UrUr_right, UrUz_right];
        b(length(P) + n_magma+1:n_magma+1:(length(P) + n_magma*(n_magma+1))) = -Dt/Ft*u2(end-n_magma+1:end)' + Bt/Ft*u1(end-n_magma+1:end)' + 1/Ft/Re.*UrR_right;

     case 'Steady'
        M(length(P) + n_magma+1:n_magma+1:(length(P) + n_magma*(n_magma+1)),:) = -[UrP_right, UrUr_right, UrUz_right];
        b(length(P) + n_magma+1:n_magma+1:(length(P) + n_magma*(n_magma+1))) = UrR_right;
end

switch BC{2}
    case 'Confined'
        j = [n_magma+1, find(z_urr>=(R/L - 1e-5))];
        M(length(P) + j,:) = 0;
        M(length(P) + j,length(P) + j) = eye(length(j));
        b(length(P) + j) = 0;

        switch BC{3}
            case 'No slip'
                if any(rem(j,n_magma+1))>0
                    'STOP'
                end
                k = [n_magma floor(j/(n_magma+1))*n_magma+n_magma + rem(j,n_magma+1)];

                M(length(P)+length(z_urr)+k,:) = 0;
                M(length(P)+length(z_urr)+k, ...
                  length(P)+length(z_urr)+k) = eye(length(k));
                b(length(P)+length(z_urr)+k) = 0;

                % Bottom right corner
                M(n_magma,:) = 0; 
                M(n_magma,n_magma) = -1;
                M(n_magma,n_magma-1) = 1/2;
                M(n_magma,2*n_magma) = 1/2;
                b(n_magma) = 0;

        end
end
Vi = V;
V = M\b;

P = V(1:length(P1))';
u = V(length(P1)+1:end)';

n = n+1;

if ~isreal(u) || ~isreal(P)
    'STOP'
end

end

P = P/(L/Mu/U);
u = u*U;

function [h1,h2,A,B,C,D,E,F] = FDcoeff(z)
h1 = max([z(2)-z(1), z(2:end-1)-z(1:end-2), z(end-1)-z(end-2)],1e-16);
h2 = max([z(3)-z(2), z(3:end)-z(2:end-1), z(end)-z(end-1)],1e-16);
A = (2*h1 + h2)./h1./(h1+h2);
B = (h1+h2)./h1./h2;
C = h1./(h1+h2)./h2;
D = h2./h1./(h1+h2);
E = (h1-h2)./h1./h2; 
F = (h1 + 2*h2)./h2./(h1+h2);