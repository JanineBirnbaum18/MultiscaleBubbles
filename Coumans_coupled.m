function [u, phi, rho, drhodt, melt_visc, eta, P, H2O, xH2O, R, Nb, T, ...
         Cc, pb_loss, m_loss, zz_p, zz_u, zz_t, t,...
         hoop_stress, along_strain_rate, transverse_strain_rate, bubble_strain_rate] = Coumans_coupled(Composition, H2Ot_0, Geometry, radius, z_int0, ...
    BC, BC_type, BC_T, flux, SolModel, DiffModel, ViscModel, EOSModel, rhoModel,  PermModel, OutgasModel, ...
    SurfTens, melt_rho,...
    rock_rho, env_rho, alpha, melt_beta, etar, Nb_0, R_0, phi_0,...
    P_0, P_f, dPdt, T_0, T_f, PTtModel, Buoyancy, dTdt, t_quench, tf,...
    solve_T, t_min, t_max, nt, n_magma)

%Gets the matlab filename
mfilename;
%Gets the folder (local computer) where the script is saved
mpath = strrep(which(mfilename),[mfilename '.m'],'');
mFolder = mpath;
%Adds the folder to the path so that local functions can be called
addpath(mFolder)

warning('off','MATLAB:illConditionedMatrix')
warning('off','MATLAB:nearlySingularMatrix')

% Gathers user-defined functions
[SolFun, DiffFun, ViscFun, m0_fun, pb_fun,...
    PTt_fun_set] = getFunctions_v2(SolModel,DiffModel, ViscModel,...
    EOSModel, PTtModel);
[SolFun, DiffFun, ViscFun, m0_fun, pb_fun,...
    PTt_fun] = getFunctions_v2(SolModel,DiffModel, ViscModel,...
    EOSModel, 'Evolving'); % 
[DynFun] = getFunctions_dynamic_dimensionless(Geometry,BC);
[DarcyFun,PermFun,WaterViscModel,...
    OutgasFun] = getFunctions_outgas(Geometry,PermModel,OutgasModel);

if solve_T
[ThermFun,cpmeltFun,cpFun,kmeltFun,kFun,...
    rhoFun] = getFunctions_thermal(Geometry,BC_type,...
    'BagdassarovDingwell1994','Bruggeman',rhoModel,alpha,T_0);
kmelt = kmeltFun(Composition);
end

%Finite difference parameters
Nodes = 50; %Number of nodes in spatial discretization

%Numerical tolerance:
%[Absolute tolerance, relative tolerance], see:
Numerical_Tolerance = [1e-5, 1e-4];
eta_max = 1e12;
Tgfun = @(T) ViscFun(H2Ot_0,T,Composition) - 1e12;
Tg = fzero(Tgfun,600+273);

%Time discretization
t = zeros([1,nt]);
dt = t_min*10; 
t(2) = dt;

% Spatial discretization
switch Geometry
    case 'Radial'
        z_t = radius*(1-logspace(0,-0.8,2*n_magma+1))*(1/(1-10^-0.8));
        z_u = z_t(1:2:end);
        z_p = z_t(2:2:end);
        g = 0;
        
    case 'Cylindrical'
        z_t = z_int0*(1-logspace(0,-0.8,2*n_magma+1))*(1/(1-10^-0.8)); 
        z_u = z_t(1:2:end);
        z_p = z_t(2:2:end);
        g = 0;%9.81;
end

% Intialize variables
[zz_t,tt_t] = meshgrid(z_t,t);
[zz_p,tt_p] = meshgrid(z_p,t);
[zz_u,tt_u] = meshgrid(z_u,t);

u = zeros(size(tt_u));
phi = ones(size(tt_p));
rho = zeros(size(tt_t));
drhodt = zeros(size(tt_p));
melt_visc = zeros(size(tt_t));
eta = zeros(size(tt_p));
beta = zeros(size(tt_p));
dbetadt = zeros(size(tt_p));
P = zeros(size(tt_p));
H2O = zeros([size(tt_p),Nodes]);
xH2O = zeros([size(tt_p),Nodes]);
R = zeros(size(tt_p));
Nb = zeros(size(tt_p));
T = zeros(size(tt_t));
Cc = zeros(size(tt_p));
pb = zeros(size(tt_p));
pb_loss = zeros(size(tt_p));
m_bub = zeros(size(tt_p));
m_loss = zeros(size(tt_p));
mean_H2O =  zeros(size(tt_t));
hoop_stress = zeros(size(tt_p));
along_strain_rate = zeros(size(tt_u));
transverse_strain_rate = zeros(size(tt_u));
bubble_strain_rate = zeros(size(tt_p));
dudr = zeros(size(tt_p));

u_t = 0*z_t;

W = Mass_SingleOxygen(Composition);
pp = 0.2*101.3e3; %2.3e3; % Partial pressure of water in surroundings

% Model tolerances
erri = 5e-4;
tol = erri;
mm = 0;
w = 0.7; % initial relaxation parameter (between 0 and 1)

% Plot options
cmap = hot(2*length(t)); % colormap

% Choose units for time axis according to total time
if tf>5*24*60*60
    t_unit = 'days';
    t_stretch = 24*60*60;
elseif tf>5*60*60
    t_unit = 'hr';
    t_stretch = 60*60;
elseif tf>5*60
    t_unit = 'min';
    t_stretch = 60;
else
    t_unit = 's';
    t_stretch = 1;
end

% Initialize plots
f3 = figure(3); clf;
switch Geometry
    case 'Radial'
        f3.Position = [400,100,600,700];
        ax3 = create_axes(5,'vertical');
    case 'Cylindrical'
        f3.Position = [400,100,1200,700];
        ax3 = create_axes(5,'horizontal');
end
set(gcf,'color','w');

labels = {'Pressure (Pa)','Velocity (m/s)','Mean H2O (wt %)','\phi','Temperature (^oC)'};
linewidth = 2;

% Initial conditions
i = 1;

while t(max([1,i-1]))<tf && i<=nt
    % Step 1
    if i == 1

        % get initial density
        phi_interp = phi_0;
        if solve_T       
            rho(i,:) = rhoFun(melt_rho,T_0 + T(1,:)).*(1-phi_interp) + density(P_0+(2.*(SurfTens)./R_0),T_0,coefficients()).*phi_interp;
        else
            rho(i,:) = melt_rho * (1-phi_interp)  + density(P_0+(2.*(SurfTens)./R_0),T_0,coefficients())*phi_interp;
        end
        
        dz = [(zz_p(1,2:end)-zz_p(1,1:end-1)),zz_p(1,end)-zz_p(1,end-1)];

        % set initial pressure
        switch Geometry
            case 'Radial'
                P(i,:) = P_0;
                Plith = P_0*ones(size(P(i,:)));
            case 'Cylindrical'

                switch Buoyancy
                    case 'True'
                        P(i,:) = cumsum(rho(i,2:2:end).*dz*g,'reverse')-rho(i,end).*dz(end)/2*g + P_0;
                        Plith = cumsum(rock_rho(i,2:2:end).*dz*g,'reverse')-rock_rho(i,end).*dz(end)/2*g + P_0;
                    case 'False'
                        P(i,:) = cumsum(rho(i,2:2:end).*dz*g,'reverse')-rho(i,end).*dz(end)/2*g + P_0;
                        Plith = cumsum(rho(i,2:2:end).*dz*g,'reverse')-rho(i,end).*dz(end)/2*g + P_0;
                end
        end

        % read other variables from inputs
        T(i,:) = T_0;

        R(i,:) = R_0;
        Nb(i,:) = Nb_0;
        phi(i,:) = phi_0;
        pb(i,:) = P(i,:)+(2.*(SurfTens)./R(i,:));
        pb_loss(i,:) = P(i,:)+(2.*(SurfTens)./R(i,:));
        m_bub(i,:) = m0_fun(R_0, pb(i,1), T(i,1));
        m_loss(i,:) = m0_fun(R_0, pb(i,1), T(i,1));

        H2O(i,:,:) = H2Ot_0;
        mean_H2O(i,:) = H2Ot_0;
        
        L = R_0*(phi_0^(-1/3) - 1); 
        xB=[0 logspace(-2,0,Nodes)*L]' + R_0;
        xH2O(i,:,:) = repmat((xB(2:end,1)+xB(1:end-1,1))/2,[1,length(z_p)])';

        melt_visc(i,:) = min(eta_max,ViscFun(H2Ot_0,T(i,:)',Composition));
        melt_visc(i,T(i)<=500) = eta_max./etar;
        eta(i,:) =  melt_visc(i,2:2:end).*etar;
        beta(i,:) = phi(i,:)./P(i,:) + (1-phi(i,:)).*melt_beta;
        dudr(i,:) = 0;
        zz_u(i,:) = z_u;
        zz_p(i,:) = z_p;
        zz_t(i,:) = z_t;

    else 

        % Adaptive time stepping
        % Allow time step to grow or shrink based on acceleration
        if i>3
            dt = dt*max(min(1.01,mean([(t(i-1)-t(i-2))/(t(i-2)-t(i-3)); (norm(u(i-2,:))/norm(u(i-1,:)))])),0.99);
        end

        % Restrict to keep spatial discretization stricly increasing
        u_t(1:2:end) = u(i-1,:);
        u_interp = interp1(zz_u(i-1,:),u(i-1,:),zz_p(i-1,:),'linear');
        u_t(2:2:end) = u_interp;
        dt_stable = max([t_max,-min([0,min(1/2*(zz_t(i-1,2:end)-zz_t(i-1,1:end-1))./(u_t(2:end)-u_t(1:end-1)))])]);

        dt = min([t_max,max([t_min,dt]),dt_stable])
        if i > 2
            if (t(i-1) - t(i-2)) < t_min
                dt = (t(i-1) - t(i-2));
            end
        end

        % If pressure is non-physical back up with a smaller time step
        if any(P(i-1,:) < 0) | any(pb_loss(i-1,:)<0) | any(isnan(P(i-1,:))) | any(abs(P(max(i-2,1),:)./P(i-1,:))>100)
            'Low pressure'
            i = i-2
            if i>2
                dt = 0.8*(t(i)-t(i-1))
            else
                i = 2
                dt = 0.8*(t(i)-t(i-1))
            end
        end

        % If spatial discretization not strictly increasing back up with a
        % smaller time step
        if any((zz_t(i-1,2:end) - zz_t(i-1,1:end-1))<0)
            'Time step too large'
            i = i-2
            if i>2
                dt = 0.8*(t(i)-t(i-1))
            else
                i = 2
                dt = 0.8*(t(i)-t(i-1))
            end
        end

        % If solution failed to converge back up with a smaller time step
        if mm>5
            'Failed to converge'
            i = i-2
            if i>2
                dt = 0.8*(t(i)-t(i-1))
            else
                i = 2
                dt = 0.8*(t(i)-t(i-1))
            end
        end

        % Caluclate coefficients for BDF2 time stepping scheme
        t(i) = t(i-1) + dt;
        if i>2
            dt2 = t(i-1) - t(i-2);
            Bt = (dt + dt2)/dt/dt2;
            Dt = dt/dt2/(dt+dt2);
            Ft = (dt2+2*dt)/dt/(dt2+dt);

        else % Bootstrap in with BDF1/backward Euler
            Dt = [0,0];
            Bt = [-1/dt,-1/dt];
            Ft = [1/dt, 1/dt];
        end

        % Main loop
        n = 0; % Is the first iteration
        m = 0; % Iteration number
        mm = 0; % Count failed tries
        P_guess = zeros(31,length(z_p));
        P_err = zeros(31,length(z_p));
        while ((erri>tol) || (n==0) || (m<3)) && (mm<5)
            % reset if stability not met
            if m>30
                dt = 0.8*dt % smaller time step
                w = max([0.2,0.9*w]) % more relaxation
                n = 0; 
                m = 0; 
                t(i) = t(i-1) + dt;
                mm = mm + 1; % Failed tries
                if i>2 % Reset time stepping coefficients after time change
                    %[~,~,~,Bt,~,Dt,~,Ft] = FDcoeff(t(i-2:i));
                    dt2 = t(i-1) - t(i-2);
                    Bt = (dt + dt2)/dt/dt2;
                    Dt = dt/dt2/(dt+dt2);
                    Ft = (dt2+2*dt)/dt/(dt2+dt);
                    else
                    Dt = [0,0];
                    Bt = [-1/dt,-1/dt];
                    Ft = [1/dt, 1/dt];
                end
            end

            % check that time step is small enough
            u_t(1:2:end) = u(i-1+n,:);
            u_t(2:2:end) = u_interp;
            dt_stable = min(1/2*(zz_t(i-1,2:end)-zz_t(i-1,1:end-1))./(u_t(2:end)-u_t(1:end-1)));
            if dt_stable<0
                dt = min([dt,-dt_stable]);
                t(i) = t(i-1) + dt;
                if i>2
                    %[~,~,~,Bt,~,Dt,~,Ft] = FDcoeff(t(i-2:i));
                    dt2 = t(i-1) - t(i-2);
                    Bt = (dt + dt2)/dt/dt2;
                    Dt = dt/dt2/(dt+dt2);
                    Ft = (dt2+2*dt)/dt/(dt2+dt);
                else
                    Dt = [0,0];
                    Bt = [-1/dt,-1/dt];
                    Ft = [1/dt, 1/dt];
                end
            end            

            % Interpolation between grids
            phi_interp = griddedInterpolant((zz_p(i-1,:)),phi(i-1+n,:),'linear','nearest');
            phi_interp = phi_interp((zz_t(i-1,:)));
            phi_interp(phi_interp>0.999) = 0.999;
            pb_interp = griddedInterpolant((zz_p(i-1,:)),pb_loss(i-1+n,:),'linear','nearest');
            pb_interp = pb_interp((zz_t(i-1,:)));

            % Bulk density
            if solve_T % Account for thermal expantion
                rho(i,:) = rhoFun(melt_rho,T(i-1+n,:)).*(1-phi_interp) + density(pb_interp,T(i-1+n,:)',coefficients()).*phi_interp;
            else
                rho(i,:) = melt_rho.*(1-phi_interp) + density(pb_interp,T(i-1+n,:)',coefficients()).*phi_interp;
            end
    
            % Thermal diffusion
            if solve_T
                cpmelt = cpmeltFun(Composition,T(i-1+n,:),mean_H2O(i-1+n,:));
                cp = cpFun(phi_interp,cpmelt,rhoFun(melt_rho,T(i-1+n,:)),T(i-1+n,:),pb_interp);
                k = kFun(phi_interp,kmelt);

                % Initial temperature solution
                PT = PTt_fun_set(P_0, P_f, dPdt,T_0,T_f,dTdt,t_quench,t(i));
  
                switch BC_type
                    case 'Forced' % Set temperature at edge
                        BC_Ti = BC_T;
    
                    case 'Dirichlet' % Set flux at edge
                        BC_Ti = PT(2);
                end
   
                if i == 2
                    T(i,:) = ThermFun(T(i-1,:),T(i-1,:),rho(i,:),cp,k,zz_t(i-1,:),dt,dt,BC_Ti,flux,'BDF1'); % Bootstrap with BDF1
                else
                    T(i,:) = ThermFun(T(i-1,:),T(i-2,:),rho(i,:),cp,k,zz_t(i-1,:),dt,t(i-1)-t(i-2),BC_Ti,flux,'BDF2'); % BDF2
                end
    
            else
                T(i,:) = PT(nt+i);
            end

            % Diffusive gas loss
            H2O_temp = squeeze(H2O(i-1,:,:));
            switch OutgasModel
                case 'Diffusive'
                    D = DiffFun([H2O(i-1,:,end),H2O(i-1,end,end)],[T(i,2:2:end),T(i,end)], [P(i-1,:),P_0], W);
                    mean_H2O_diff = OutgasFun([H2O(i-1,:,end),SolFun(BC_T(end),pp)],[H2O(i-1,:,end),SolFun(BC_T(end),pp)],D,[zz_p(i-1,:),zz_u(i-1,end)],dt,dt,SolFun(BC_T(end),pp),'BDF1');
                    H2O_temp(:,end) = mean_H2O_diff(1:end-1);     
            end

            % Run Coumans 2020 for each pressure node
            for j = 1:length(z_p)

                % Skip nodes that can't grow
                if (R(i-1,j)<=1.01e-6 && SolFun(T(i-1,2*j),pb_loss(i-1,j))>mean_H2O(i-1,2*j)) || T(i-1,2*j)<Tg
                    Nb(i,j) = Nb(i-1,j);
                    R(i,j) = R(i-1,j);
                    phi(i,j) = phi(i-1,j);
                    H2O(i,j,:) = H2O(i-1,j,:);

                    switch OutgasModel
                        case 'Diffusive'
                            D = DiffFun(H2O(i-1,j,:),T(i,j) + 0*xH2O(i,j,:), P(i-1,j) + 0*xH2O(i,j,:), W);
                            H2O(i,j,:) = OutgasFun([squeeze(H2O(i-1,j,1:end-1))',SolFun(BC_T(end),pp)],[squeeze(H2O(i-1,j,1:end-1))',SolFun(BC_T(end),pp)],squeeze(D)',squeeze(xH2O(i-1,j,:))',dt,dt,SolFun(BC_T(end),pp),'BDF1');
                    end
                    
                    xH2O(i,j,:) = xH2O(i-1,j,:);
                    mean_H2O(i,2*j) = trapz(squeeze(xH2O(i,j,:)).^3,squeeze(H2O(i,j,:)))./(xH2O(i,j,end).^3-xH2O(i,j,1).^3);
                    pb(i,j) = pb(i-1,j);
                    pb_loss(i,j) = pb_loss(i-1,j);
                    m_bub(i,j) = m_bub(i-1,j);
                    m_loss(i,j) = m_loss(i-1,j);

                else         
                    % Project pressure and temperature
                    dTdti = (T(i,j) - T(i-1,j))/dt;
                    if i > 2
                        if n == 0 % If first iteration at timestep, use current slope
                            if i>3
                                h2 = t(i-1) - t(i-2);
                                h1 = t(i-2) - t(i-3);
                                D = h2./h1./(h1+h2);
                                B = (h1+h2)./h1./h2;
                                F = (h1 + 2*h2)./h2./(h1+h2);
                                dPdti = D*P(i-3,j)-B*P(i-2,j)+F*P(i-1,j) + ...
                                    (2*h2/(h1*h2*(h1+h2))*P(i-3,j)-2*(h1+h2)./(h1.*h2.*(h1+h2))*P(i-2,j) + ...
                                    2*h1/(h1*h2*(h1+h2))*P(i-1,j))*dt/2;
                            else
                                dPdti = 0; % Start from rest
                            end
                            Pf = max(1e-3,P(i-1,j) + dPdti*dt);
                        else % Otherwise use most recent guess
                            Pf = P(i,j);
                            dPdti = (P(i,j)-P(i-1,j))/dt;
                        end
                    else % Bootstrap in with a diffusion-based projection
                        if n == 0
                            H2O_forecast = (SolFun(T_0,P_0)-H2Ot_0).*erfc(xH2O(1,:,:)./(2*sqrt(0.1*DiffFun(H2Ot_0,T_0,P_0,W)*(t(i))))) + H2Ot_0;
                            I1=trapz(squeeze(xH2O(1,j,:)),squeeze((1/100)*H2O(1,j,:)).*squeeze(xH2O(1,j,:).^2),1);
                            I2=trapz(squeeze(xH2O(1,j,:)),squeeze((1/100)*H2O_forecast(1,j,:)).*squeeze(xH2O(1,j,:).^2),1);            
                            m_forecast=m_loss(1,j)+4.*pi.*melt_rho.*(I1-I2);

                            Pf = w*pb_fun(m_forecast,T_0,R_0) + (1-w)*P_0;
                            dPdti = (Pf - P(i-1,j))./dt;
                        else
                            dPdti = (P(i,j) - P(i-1,j))./dt;
                            Pf = P(i,j);
                        end
                    end

                    %if (rem(m,10)==0) & m>1
                    %    fitobject = fit(P_err(2:m,j),P_guess(2:m,j),'linear','Exclude',P_guess(2:m,j)<1);
                    %    Pf = fitobject(0);
                    %    dPdti = (Pf - P(i-1,j))./dt;
                    %end

                    P_guess(m+1,j) = Pf;
                    
                    % Call Coumans 2020
                    [ti, Ri, phii, Pii, Tii, x_out, H2Ot_all, Nbi, pbi, mi] =  Numerical_Model_v2(Composition, SolModel, DiffModel,...
                        ViscModel, EOSModel,OutgasModel, 'Evolving', SurfTens, melt_rho, Nodes,...
                        R(i-1,j),...
                        [squeeze(H2O(1,j,:)),H2O_temp(j,:)', squeeze(xH2O(i-1,j,:))], m_loss(i-1,j),...
                        Nb(i-1,j), 0, dt, T(i-1,j), T(i,j), dTdti,...
                        P(i-1,j), Pf,...
                        dPdti, 1e10, Numerical_Tolerance,...
                        eta(i-1+n,:), zz_p(i-1,:), j, Geometry, radius);

                    % Update solution with SOR

                    Nb(i,j) = ((1-w)*(1-n) + w)*Nbi(end) + n*(1-w)*Nb(i,j);
                    R(i,j) = ((1-w)*(1-n) + w)*Ri(end) + n*(1-w)*R(i,j);
                    phi(i,j) = ((1-w)*(1-n) + w)*phii(end) + n*(1-w)*phi(i,j);
                    xH2O(i,j,:) = ((1-w)*(1-n) + w)*squeeze(x_out(:,end)) + n*squeeze((1-w)*xH2O(i,j,:));
                    H2O(i,j,:) = ((1-w)*(1-n) + w)*squeeze(H2Ot_all(:,end)) + n*squeeze((1-w)*H2O(i,j,:));
                    mean_H2O(i,2*j) = trapz(squeeze(xH2O(i,j,:)).^3,squeeze(H2O(i,j,:)))./(xH2O(i,j,end).^3-xH2O(i,j,1).^3);
                    pb(i,j) = ((1-w)*(1-n) + w)*pbi(end) + n*(1-w)*pb(i,j);
                    m_bub(i,j) = ((1-w)*(1-n) + w)*mi(end) + n*(1-w)*m_bub(i,j);
                end
            end

            % Interpolate water concentration to all nodes for viscosity
            % estimate
            switch OutgasModel
                case 'Diffusive'
                    H2O_diff = griddedInterpolant([zz_p(i-1,:), zz_u(i-1,end)],[mean_H2O(i,2:2:end),SolFun(BC_T(end),pp)],'linear','nearest');
                    mean_H2O(i,1:2:end) = H2O_diff([zz_u(i-1,1:end-1),(zz_p(i-1,end)+zz_u(i-1,end))/2]);
                case 'None'
                    H2O_diff = griddedInterpolant((zz_p(i-1,:)),mean_H2O(i,2:2:end),'linear','nearest');
                    mean_H2O(i,1:2:end) = H2O_diff((zz_u(i-1,:)));
            end
            
            % Interpolate between grids
            phi_interp = griddedInterpolant((zz_p(i-1,:)),phi(i,:),'linear','linear');
            phi_interp = phi_interp((zz_t(i-1,:)));
            phi_interp(phi_interp>0.999) = 0.999;
            phi_interp(phi_interp<phi_0) = phi_0; 
            pb(i,pb(i,:)<1) = 1;
            pb_interp = griddedInterpolant((zz_p(i-1,:)),pb(i,:),'linear','linear');
            pb_interp = pb_interp((zz_t(i-1,:)));
            pb_interp(pb_interp<1) = 1;
            
            % Update bulk material properties
            gas_rho = density(pb_interp,T(i,:)',coefficients());
            if solve_T
                rho(i,:) = rhoFun(melt_rho,T(i,:)).*(1-phi_interp) + gas_rho.*(phi_interp);
            else
                rho(i,:) = melt_rho*(1-phi_interp) + gas_rho.*(phi_interp);
            end

            melt_visc(i,:) = min(eta_max,ViscFun(mean_H2O(i,:)',T(i,:)',Composition));
            melt_visc(i,T(i,:)<=500) = eta_max./etar;
            beta(i,:) = phi(i,:)./pb(i,:) + (1-phi(i,:)).*melt_beta;

            % And their time derivatives
            if i == 2
                drhodt(i,:) = (rho(i,2:2:end) - rho(i-1,2:2:end))./dt;
                dbetadt(i,:) = (beta(i,:) - beta(i-1,:))./dt;
            else
                drhodt(i,:) = Dt*rho(i-2,2:2:end) - Bt*rho(i-1,2:2:end) + Ft*rho(i,2:2:end);
                dbetadt(i,:) = Dt*beta(i-2,:) - Bt*beta(i-1,:) + Ft*beta(i,:);
            end
            

            % Viscous flow
            if i == 2 
                [Pi, ui] = DynFun(P(i-1,:)-mean(Plith),u(i-1,:), ...
                    P(i-1,:)-mean(Plith),u(i-1,:),Plith-mean(Plith), ...
                    rho(i,:),drhodt(i,:),phi(i,:),R(i,:),SurfTens,melt_visc(i,:).*etar, ...
                    beta(i,:),dbetadt(i,:),radius,g,zz_p(i-1,:),zz_u(i-1,:),dt,dt,'BDF1');
            else
                [Pi, ui] = DynFun(P(i-1,:)-mean(Plith),u(i-1,:), ...
                    P(i-2,:)-mean(Plith),u(i-2,:),Plith-mean(Plith), ...
                    rho(i,:),drhodt(i,:),phi(i,:),R(i,:),SurfTens,melt_visc(i,:).*etar, ...
                    beta(i,:),dbetadt(i,:),radius,g,zz_p(i-1,:),zz_u(i-1,:),dt,t(i-1)-t(i-2),'BDF2');
            end

            P_err(m+1,:) = P(i,j) - (Pi + mean(Plith));
            P(i,:) = n*(1-w)*P(i,:) + ((1-n)*(1-w) + w)*(Pi + mean(Plith));

            erri = norm(u(i,:) - ui)./max(((norm(u(i,:)) + norm(u(i-1,:)))/2),1e-12);
            u(i,:) = n*(1-w)*u(i,:) + ((1-n)*(1-w) + w)*ui;
            u_interp = interp1(zz_u(i-1,:),u(i,:),zz_p(i-1,:),'linear');

            % Calculate capillary number
            switch Geometry
                case 'Radial'
                    dudr(i,:) = 0*(u(i,2:end) - u(i,1:end-1))./(zz_u(i-1,2:end)-zz_u(i-1,1:end-1));
                case 'Cylindrical'
                    dudr(i,:) = 3*u_interp/radius;
            end

            if i>2
                d2udrdt = (Ft.*dudr(i,:) - Bt.*dudr(i-1,:) + Dt*dudr(i-2,:));
            else
                d2udrdt = 0;
            end

            Cc(i,:) = max(sqrt((d2udrdt./dudr(i,:)).^2 + dudr(i,:).^2),1e-10).*R(i,:).*melt_visc(i,2:2:end).*etar/SurfTens;
            eta0 = (1-phi(i,:)).^(-1);
            etainf = (1-phi(i,:)).^(5/3);
            eta(i,:) = melt_visc(i,2:2:end).*etar.*(etainf + (eta0-etainf)./(1+(6/5*Cc(i,:)).^2));
            eta(eta>eta_max) = eta_max;

            % Permeable outgassing
            
            switch PermModel
                case 'None'
                    m_loss(i,:) = m_bub(i,:); 
                    pb_loss(i,:) = pb(i,:);
                otherwise
                    gas_rho = density(Plith(end),T(i,2:2:end)',coefficients());
                    m_loss(i,:) = n*(1-w)*m_loss(i,:) + ((1-n)*(1-w) + w)*DarcyFun(m0_fun,m_bub(i,:),pb(i,:),Plith,radius,z_p,z_u,...
                    PermFun(phi(i,:),Cc(i,:)),...
                    WaterViscModel(gas_rho,T(i,2:2:end)),...
                    gas_rho,Nb(i,:),R(i,:),T(i,2:2:end),dt,gas_rho);
            
                    pb_loss(i,:) = pb_fun(m_loss(i,:), T(i,2:2:end), R(i,:));
            end

            % Calculate failure
            hoop_stress(i,:) = (P(i,:)-P_0).*(radius)./(2*((zz_u(i-1,end)-zz_p(i-1,:))));
            along_strain_rate(i,2:end) = (u(i,2:end) - u(i,1:end-1))./(zz_u(i-1,2:end)-zz_u(i-1,1:end-1)).*melt_visc(i,3:2:end)/1e10;
            switch Geometry
                case 'Radial'
                    transverse_strain_rate(i,2:end) = 1./zz_u(i-1,2:end).*u(i,2:end).*melt_visc(i,3:2:end)/1e10;
                case 'Cylindrical'
                    transverse_strain_rate(i,:) = (3*u(i,:)/radius).*melt_visc(i,1:2:end)/1e10;
            end
            
            if i>2
                dRdt = (Ft.*R(i,:) - Bt.*R(i-1,:) + Dt*R(i-2,:));
            else
                dRdt = (R(i,:)-R(i-1,:))./dt;
            end
            bubble_strain_rate(i,:) = 1./xH2O(i,:,1).*(abs(dRdt)).*ViscFun(squeeze(H2O(i,:,1)),T(i,2:2:end),Composition)./1e10;

            % Update coordinates
            if i == 2
                zz_p(i,:) = zz_p(i-1,:)+u_interp.*dt;
                zz_u(i,:) = zz_u(i-1,:)+u(i,:).*dt;
    
                zz_t(i,1:2:end) = zz_t(i-1,1:2:end)+u(i,:).*dt;
                zz_t(i,2:2:end) = zz_t(i-1,2:2:end)+u_interp.*dt;
    
            else
                zz_p(i,:) = -Dt/Ft*zz_p(i-2,:) + Bt/Ft*zz_p(i-1,:) + 1/Ft.*u_interp;
                zz_u(i,:) = -Dt/Ft*zz_u(i-2,:) + Bt/Ft*zz_u(i-1,:) + 1/Ft.*u(i,:);
    
                zz_t(i,1:2:end) = -Dt/Ft*zz_t(i-2,1:2:end) + Bt/Ft*zz_t(i-1,1:2:end) + 1/Ft.*u(i,:);
                zz_t(i,2:2:end) = -Dt/Ft*zz_t(i-2,2:2:end) + Bt/Ft*zz_t(i-1,2:2:end) + 1/Ft.*u_interp;
            end
    
            switch Geometry
                case 'Radial'
                    radius = zz_t(i-1,end);
            end

            % Plot results
            figure(8);
            if m == 0
                clf;
            end
            subplot(5,2,1); hold on;
            plot(m,mean(P(i,:)),'o','MarkerFaceColor','k')
            ylabel('Pressure (Pa)')
            title(['i = ' string(i)])
            subplot(5,2,2); hold on;
            plot(m,mean(drhodt(i,:)),'o','MarkerFaceColor','k')
            title(['t = ' string(t(i))])
            ylabel('d\rho/dt (m/s)')
            subplot(5,2,3); hold on;
            plot(m,mean(phi_interp),'o','MarkerFaceColor','k')
            ylabel('Vesicularity')

            subplot(5,2,4); hold on;
            plot(m,mean(mean_H2O(i,:)),'o','MarkerFaceColor','k')
            ylabel('H2O')
            subplot(5,2,5); hold on;
            plot(m,mean(pb_interp),'o','MarkerFaceColor','k')
            ylabel('Bubble pressure')
            subplot(5,2,6); hold on;
            plot(m,mean(beta(i,:)),'o','MarkerFaceColor','k')
            ylabel('Compressibility')

            subplot(5,2,7); hold on;
            plot(m,mean(eta(i,:)),'o','MarkerFaceColor','k')
            ylabel('Viscosity')
            subplot(5,2,8); hold on;
            plot(m,mean(u(i,:)),'o','MarkerFaceColor','k')
            ylabel('Velocity')

            subplot(5,2,9); hold on;
            plot(m,mean(Nb(i,:)),'o','MarkerFaceColor','k')
            ylabel('Nb')
            subplot(5,2,10); hold on;
            plot(m,mean(R(i,:)),'o','MarkerFaceColor','k')
            ylabel('Radius')

            n = 1;
            m = m+1;
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%% figure 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(3)
    
    if rem(i-1,round((nt-1)/10))==0
        ys = {P-P_0,u,mean_H2O,phi,T-273.15};
        xs = {zz_p,zz_u,zz_t,zz_p,zz_t};
        switch Geometry
            case 'Radial'
                for j = 1:5
                    set(f3,'CurrentAxes',ax3(j))
                    hold on;
                    h3 = plot(xs{j}(i,:),ys{j}(i,:),'color',cmap(i,:),'linewidth',linewidth);
                    xlim([0,max(zz_t(:))*1.01])
        
                    ylabel(labels(j))
                    box('on')
                end
                xlabel('Radius (m)')
                for j=1:4
                    set(f3,'CurrentAxes',ax3(j))
                end

            case 'Cylindrical'
                for j = 1:5
                     set(f3,'CurrentAxes',ax3(j))
                     hold on;
                     h3 = plot(ys{j}(i,:),xs{j}(i,:),'color',cmap(i,:),'linewidth',linewidth);
                     ylim([0,max(zz_t(:))*1.01])
        
                    xlabel(labels(j))
                    box('on')
                end
                set(f3,'CurrentAxes',ax3(1))
                ylabel('Height (m)')
        
                for j=2:5
                    set(f3,'CurrentAxes',ax3(j))
                    %yticklabels([])
                end
        end
    drawnow
    end
    
    i = i + 1
end

if t(i-1)>tf
    u = u(1:i-1,:);
    phi = phi(1:i-1,:);
    rho = rho(1:i-1,:);
    drhodt = drhodt(1:i-1,:);
    melt_visc = melt_visc(1:i-1,:);
    eta = eta(1:i-1,:);
    P = P(1:i-1,:);
    H2O = H2O(1:i-1,:,:);
    xH2O = xH2O(1:i-1,:,:);
    mean_H2O = mean_H2O(1:i-1,:);
    R = R(1:i-1,:);
    Nb = Nb(1:i-1,:);
    T = T(1:i-1,:);
    Cc = Cc(1:i-1,:);
    pb = pb(1:i-1,:);
    pb_loss = pb_loss(1:i-1,:);
    zz_p = zz_p(1:i-1,:);
    zz_u = zz_u(1:i-1,:);
    zz_t = zz_t(1:i-1,:);
    t = t(1:i-1);
end

colormap(cmap(1:end/2,:))
c = colorbar('Position',[0.93 0.168 0.015 0.7]);
c.Label.String = ['Time (' t_unit ')'];
clim([0,max(t)/t_stretch]);  

%%%%%%%%%%%%%%%%%%%%%%%% figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z2 = n_magma-1;

f4 = figure(4); clf;
f4.Position = [400,100,600,700];
set(gcf,'color','w');

ys = {P-P_0,pb_loss-P-(2.*(SurfTens)./R),mean_H2O,phi};
labels = {'\Delta Pressure (Pa)','Bubble overpressure (Pa)','Mean H2O (wt %)','\phi'};
ax4 = create_axes(4,'vertical');
for j = 1:4
    set(f4,'CurrentAxes',ax4(j))
    hold on;
    for i = 1:length(t)
        h1 = plot(t/t_stretch,ys{j}(:,2),'b','linewidth',linewidth);
        h2 = plot(t/t_stretch,ys{j}(:,z2),'r','linewidth',linewidth);
    end

    ylabel(labels(j))
    box('on')
end

xlabel(['Time (' t_unit ')'])
for j=1:3
    set(f4,'CurrentAxes',ax4(j))
    xticklabels([])
end
plot(t/t_stretch,SolFun(T(:,2),P(:,2)),'b--','linewidth',linewidth);
plot(t/t_stretch,SolFun(T(:,z2),P(:,z2)),'r--','linewidth',linewidth);

legend([h1,h2],'z=' + string(z_p(2)) + ' m',...
    'z=' + string(z_p(z2)) + ' m')

%%%%%%%%%%%%%%%%%%%%% figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%
f5 = figure(5); clf;
f5.Position = [400,100,600,700];
set(gcf,'color','w');

ys = {mean_H2O,R,Nb, phi};
labels = {'H2O (wt %)','Bubble radius (m)','Bubble number density','\phi'};

ax5 = create_axes(4,'vertical');
for j = 1:4
    set(f5,'CurrentAxes',ax5(j))
    hold on;
    for i = 1:length(t)
        h1 = plot(t/t_stretch,ys{j}(:,2),'b','linewidth',linewidth);
        h2 = plot(t/t_stretch,ys{j}(:,z2),'r','linewidth',linewidth);
    end

    ylabel(labels(j))
    box('on')
end
xlabel(['Time (' t_unit ')'])
for j=1:3
    set(f5,'CurrentAxes',ax5(j))
    xticklabels([])
end

legend([h1,h2],'z=' + string(z_p(2)) + ' m',...
    'z=' + string(z_p(z2)) + ' m')

%%%%%%%%%%%%%%%%%%%%%%%%% figure 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f6 = figure(6); clf;
f6.Position = [400,100,600,700];
set(gcf,'color','w');

ys = {phi,u,Cc};
labels = {'\phi','Velocity (m/s)','Cc'};

ax6 = create_axes(3,'vertical');
for j = 1:3
    set(f6,'CurrentAxes',ax6(j))
    hold on;
    for i = 1:length(t)
        h1 = plot(t/t_stretch,ys{j}(:,2),'b','linewidth',linewidth);
        h2 = plot(t/t_stretch,ys{j}(:,z2),'r','linewidth',linewidth);
    end

    ylabel(labels(j))
    box('on')
end
xlabel(['Time (' t_unit ')'])
for j=1:2
    set(f6,'CurrentAxes',ax6(j))
    xticklabels([])
end

set(ax6(3),'YScale','log')

legend([h1,h2],'z=' + string(z_p(2)) + ' m',...
    'z=' + string(z_p(z2)) + ' m')

%%%%%%%%%%%%%%%%%%%%%%%%% figure 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f7 = figure(7); clf;
f7.Position = [400,100,600,700];
set(gcf,'color','w');

ys = {P-P_0,u,T-273.15,phi};
labels = {'\Delta P (Pa)','Velocity (m/s)','Temperature (^oC)','Vesciularity'};

ax7 = create_axes(4,'vertical');
for j = 1:4
    set(f7,'CurrentAxes',ax7(j))
    hold on;
    for i = 1:length(t)
        h1 = plot(t/t_stretch,ys{j}(:,2),'b','linewidth',linewidth);
        h2 = plot(t/t_stretch,ys{j}(:,z2),'r','linewidth',linewidth);
    end

    ylabel(labels(j))
    box('on')
end
xlabel(['Time (' t_unit ')'])
for j=1:3
    set(f7,'CurrentAxes',ax7(j))
    xticklabels([])
end

legend([h1,h2],'z=' + string(z_p(2)) + ' m',...
    'z=' + string(z_p(z2)) + ' m')
end

% Supporting functions

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
end

function rho = density(P,T,b)
% convert P to bars from input (pascals)
P = P.*1e-5; %1 pascal a 1e-5 bars

a=zeros(10,length(T));
for i=1:10
    a(i,:)=b(i,1).*T.^-4 + b(i,2).*T.^-2 + b(i,3).*T.^-1 +...
        b(i,4) + b(i,5).*T + b(i,6).*T.^2;
end
% PRT = P/RT where P [bars], R [cm^3*bar/K/mol] and T [K]
PRT = P./(83.14472*T);

rho = 0*T';
% solve implicit equation for rho and convert to kg/m^3
for i=1:length(T)
rho(i) = fzero(@PS_myfun,0.001,[],a(:,i),PRT(i))*18.01528*1000;
end
end

% the function from Pitzer & Sterner 1994, which takes the matrix of
% coefficients a and P/RT as arguments; rho is a first guess for the
% density [g/mol]
function y = PS_myfun(rho,a,PRT)
y = (rho+a(1)*rho^2-rho^2*((a(3)+2*a(4)*rho+3*a(5)*rho^2+4*a(6)*rho^3)/((a(2)+a(3)*rho+a(4)*rho^2+a(5)*rho^3+a(6)*rho^4)^2))+a(7)*rho^2*exp(-a(8)*rho)+a(9)*rho^2*exp(-a(10)*rho)) - PRT;
end

function [h1,h2,A,B,C,D,E,F] = FDcoeff(z)
h1 = [z(2)-z(1), z(2:end-1)-z(1:end-2), z(end-1)-z(end-2)];
h2 = [z(3)-z(2), z(3:end)-z(2:end-1), z(end)-z(end-1)];
A = (2*h1 + h2)./h1./(h1+h2);
B = (h1+h2)./h1./h2;
C = h1./(h1+h2)./h2;
D = h2./h1./(h1+h2);
E = (h1-h2)./h1./h2; 
F = (h1 + 2*h2)./h2./(h1+h2);
end

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
end