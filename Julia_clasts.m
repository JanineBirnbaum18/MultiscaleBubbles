close all
clear
warning('off','MATLAB:illConditionedMatrix')
warning('off','MATLAB:nearlySingularMatrix')

%Melt composition
%[SiO2 TiO2 Al2O3 FeO(T) MnO MgO CaO Na2O K2O P2O5 H2O F2O-1]
Composition = [75.17 0.22 12.02 3.13 0.11 0.09 1.66 4.58 2.88 0 0 0];
H2Ot_0 = 0.115; %initial water concentration (wt. %)
%Dynamics geometry
Geometry = 'Radial'; 
BC = 'Symmetry';
BC_type = 'Dirichlet';
BC_T = 1006+273.15;
flux = 0;

%Material properties
SolModel = 'Schunke';
DiffModel = 'Zhang 2010 Metaluminous';
ViscModel = 'Hess and Dingwell 1996';

%Model options
EOSModel = 'Pitzer and Sterner';

%Thermal properties
rhoModel = 'BagdassarovDingwell1994';
alpha = 0;

% Outgassing options
PermModel = 'None';
OutgasModel = 'Diffusive';

%Constants used for calculations
SurfTens = 0.08; %Value for surface tension (N/m) Shea 2017: 0.079
melt_rho = 2400; %Melt density in kg/m^3
rock_rho = 2400; %Country rock density in kg/m^3
env_rho = 1; %Surrounding air or water density in kg/m^3
melt_beta = 2.6e-11; % Malfait et al. (2011)

%Spatial parameters
Nb_0 = 5e10; %Bubble number density (number per m^3)
R_0 = 1e-5; %R_0 (m) set independently
phi_0 = (4/3).*pi()*R_0.^3./(1/Nb_0);

% P-T pathway 
P_0 = 101e3; % Initial surface pressure
P_f = 101e3; % Final surface pressure
dPdt = 1e-10;
T_0 = 806 + 273.15; %Initial temperature in K
T_f = 1006 + 273.15;
PTtModel = 'Jenny';
Buoyancy = 'False';
dTdt = [10/60 -10/60];
t_quench = 10*60*60+30*60;
tf = 1e5;
solve_T = true;

% Discretization
nt = 200;
n_magma = 10;
t_min = 20;
t_max = 600;

radii = [0.5, 3, 6]*1e-3;

for i = 1:length(radii)
    radius = radii(i);
    z_int0 = radii(i);
    [u, phi, rho, drhodt, melt_visc, eta, P, H2O, xH2O, R, Nb, T, Cc, pb,...
    m_loss, zz_p, zz_u, zz_t, t] = Coumans_coupled(Composition, H2Ot_0,...
         Geometry, radius, z_int0, ...
    BC, BC_type, BC_T, flux, SolModel, DiffModel, ViscModel, EOSModel, rhoModel,...
    PermModel, OutgasModel, SurfTens, melt_rho,...
    rock_rho, env_rho, alpha, melt_beta, 1, Nb_0, R_0, phi_0,...
    P_0, P_f, dPdt, T_0, T_f, PTtModel, Buoyancy, dTdt, t_quench, tf,...
    solve_T, t_min, t_max, nt, n_magma);

    figure(9); hold on; 
    plot(t+78*60,sum((zz_u(:,2:end).^3 - zz_u(:,1:end-1).^3).*phi,2)./(zz_u(:,end)).^3,'DisplayName','R = ' + string(radius*1e3) + ' mm');
end
set(gca,'XScale','log')
xlim([1e3,2e5])
ylim([0,0.8])

cmap = hot(2*length(t));

t_unit = 's';
if contains(t_unit,'s')
    t_stretch = 1;
elseif contains(t_unit,'min')
    t_stretch = 60;
elseif contains(t_unit,'hr')
    t_stretch = 60*60;
elseif contains(t_unit,'days')
    t_stretch = 60*60*24;
end

f1 = figure(1);
f1.Position = [400,100,600,700];
ax1 = create_axes(4,'vertical');
set(gcf,'color','w');

labels = {'\phi','\Delta Pressure (Pa)','Velocity (m/s)','Temperature (^oC)'};
linewidth = 2;
ys = {phi,P-101e3,u,T-273};
xs = {zz_p,zz_p,zz_u,zz_t};

for i=1:max([1, floor(length(t)/10)]):length(t)
    for j = 1:4
        set(f1,'CurrentAxes',ax1(j))
        hold on;
        h1 = plot(xs{j}(i,:),ys{j}(i,:),'color',cmap(i,:),'linewidth',linewidth);

        ylabel(labels(j))
        box('on')
    end
    xlabel('Radius (m)')
    for j=1:3
        set(f1,'CurrentAxes',ax1(j))
        xticklabels([])
    end
end


dat = struct('t',t,'R',R,'R0',R_0,'Nb',Nb,'Nb0',Nb_0,'phi',phi,'Pb',pb,'m',m_loss, ...
            'P',P,'P0',P_0,'dPdt',dPdt,'u',u, ...
            'T',T,'T0',T_0,'Tf' ,T_f,'dTdt',dTdt,'t_quench',t_quench,'PModel',PTtModel, ...
            'xH2O',xH2O ,'H2O',H2O, 'H2O0',H2Ot_0, ...
            'zz_p',zz_p ,'zz_u',zz_u, 'zz_t',zz_t, ...
            'rho',rho,'drhodt',drhodt, 'melt_visc',melt_visc,'eta',eta, 'Cc',Cc, ...
            'Composition',Composition,'DiffModel_H2O',DiffModel,...
            'EOSModel',EOSModel,'SolModel',SolModel ,...
            'ViscModel',ViscModel,'SurfTens',SurfTens,'melt_rho',melt_rho);
    save(['./Results/Radius' num2str(radius*1e-6,'%.0f') '.mat'],'dat')


