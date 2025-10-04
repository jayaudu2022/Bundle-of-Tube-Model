%% Bundle of Tubes (BOT) Salt Precipitation Model - Clean Version
clear; clc;

%% 1. Physical Parameters
params.T      = 323.15;
params.Pg     = 101325;
params.RH     = 0.0;
params.rho    = 1000;
params.sigma  = 0.0718;
params.mu     = 1.0142e-3;
params.costh  = cosd(30);
params.g      = 9.81;
params.rho_s  = 2165;
params.Vm     = 1.802e-5;
params.Rg     = 8.314;
params.Psat   = 12.33e3;
params.Mw     = 0.018;
params.Dv     = 2.5e-4;
params.L      = 3.0;

% Salt parameters
params.m0     = 6.0;
params.m_sat  = 6.14;
params.M_NaCl = 0.05844;

% Geometry
params.h_fixed = 5e-3;
params.r_min   = 5e-6;

% Precipitation
params.l_rim       = 1e-3;
params.z_min       = 1e-6;
params.use_gravity = 1;
params.output_interval = inf;  

%% 2. BOT-Specific Parameters
params.N_tubes      = 10;
r_ch                = 85e-6;
s                   = 0.3;
params.r_min_cutoff = 50e-6;
params.r_max_cutoff = 200e-6;

fprintf(' BOT SIMULATION \n');
fprintf('Tubes: %d, Range: %.0f-%.0f μm\n', ...
        params.N_tubes, params.r_min_cutoff*1e6, params.r_max_cutoff*1e6);

%% 3. Generate Tube Distribution
rng(42);
r_tubes = sample_trunc_lognormal(params.N_tubes, r_ch, s, params.r_min_cutoff, params.r_max_cutoff);
h_tubes = params.h_fixed * ones(params.N_tubes, 1);
params.V_pore_initial = sum(pi * r_tubes.^2 .* h_tubes);

%% 4. Initial Conditions
y0 = [];
for i = 1:params.N_tubes
    V0_i       = pi * r_tubes(i)^2 * h_tubes(i);
    M_water_0  = params.rho * V0_i;
    Md0_i      = params.m0 * M_water_0 * params.M_NaCl;
    y0_i       = [h_tubes(i); Md0_i; r_tubes(i); 0; 0; M_water_0];
    y0         = [y0; y0_i];
end

%% 5. Run Simulation
tspan   = [0 1e8];
options = odeset('RelTol',1e-5,'AbsTol',1e-10,'MaxStep',1,'Events',@(t,y) bot_closure_events(t,y,params));

fprintf('Running simulation...\n');
tic;
[t, y] = ode15s(@(t,y) bot_odefun(t,y,params), tspan, y0, options);
fprintf('Completed in %.1f seconds\n', toc);

%% 6. Save Results
save('results_BOT.mat','t','y','params','r_tubes');
fprintf('Saved as results_BOT.mat\n');

%%  FUNCTIONS 

function r = sample_trunc_lognormal(N_tubes, r_ch, s, rmin, rmax)
    mu = log(r_ch);
    Fa = 0.5*(1 + erf((log(rmin) - mu)/(s*sqrt(2))));
    Fb = 0.5*(1 + erf((log(rmax) - mu)/(s*sqrt(2))));
    U  = Fa + (Fb - Fa)*rand(N_tubes,1);
    r  = exp(mu + s*sqrt(2)*erfinv(2*U - 1));
    r  = min(max(r,rmin),rmax);
end

function dydt = bot_odefun(t, y, params)
    persistent step_count
    if isempty(step_count), step_count = 0; end
    step_count = step_count + 1;
    
    N_tubes = params.N_tubes;
    state_size = 6;
    dydt = zeros(size(y));
    
    for i = 1:N_tubes
        idx_start = (i-1)*state_size + 1;
        idx_end   = i*state_size;
        y_tube    = y(idx_start:idx_end);
        
        if y_tube(3) > params.r_min
            dydt_tube = single_tube_ode(t, y_tube, params);
        else
            dydt_tube = zeros(6,1);
        end
        
        dydt(idx_start:idx_end) = dydt_tube;
    end
    
    if mod(step_count, 2000) == 0
        active = sum(y(3:6:end) > params.r_min);
        fprintf('t=%.0fs: %d/%d active\n', t, active, N_tubes);
    end
end

function dydt = single_tube_ode(~, y, params)
    h = max(y(1),0);
    Md = max(y(2),0);
    r  = max(y(3), params.r_min);
    M_water = y(6);

    dydt = zeros(6,1);

    h_max = 0.9 * params.L;
    h_eff = min(h,h_max);
    A = pi*r^2;
    V = A*h_eff;

    Pc    = 2*params.sigma*params.costh/r;
    pgrav = params.use_gravity*params.rho*params.g*h_eff;
    Rh    = 8*params.mu*h_eff/(pi*r^4);

    kelv = min(2*params.sigma*params.Vm*params.costh/(r*params.Rg*params.T),50);
    f1 = params.RH*params.Psat;
    f2 = params.Psat*exp(kelv);
    y1 = min(f1/params.Pg,0.999999);
    y2 = min(f2/params.Pg,0.999999);

    if (1-y2)<=1e-10
        lnterm = 50;
    else
        lnterm = max(log((1-y1)/(1-y2)),0);
    end

    L_eff = max(params.L - h_eff, 1e-9);
    Jm = (params.Pg*params.Dv*params.Mw)/(params.Rg*params.T) * (lnterm/L_eff);
    E  = Jm*A/params.rho;

    if h>=h_max
        Q = min((Pc-pgrav)/max(Rh,1e-30), E);
        dydt(1)=0;
    else
        Q = (Pc-pgrav)/max(Rh,1e-30);
        dydt(1) = (Q-E)/A;
    end

    Md_inflow = params.m0*Q*params.rho*params.M_NaCl;
    mol_current = (Md/params.M_NaCl)/max(M_water,1e-20);

    if mol_current>params.m_sat
        supersat = mol_current-params.m_sat;
        k_precip_rate = 1e-4;
        if mol_current>1.5*params.m_sat
            k_precip_rate = k_precip_rate*100;
        elseif mol_current>1.2*params.m_sat
            k_precip_rate = k_precip_rate*10;
        end
        R_precip_mol = k_precip_rate*supersat*M_water;
        R_precip = min(R_precip_mol*params.M_NaCl, 0.9*Md);
    else
        R_precip = 0;
    end
    R_precip = max(R_precip,0);

    if R_precip>0
        z_dep = min(max(h_eff,0), max(params.z_min, params.l_rim));
        dr_dt = -R_precip/(params.rho_s*2*pi*r*z_dep);
    else
        dr_dt=0;
    end

    dydt(2) = Md_inflow - R_precip;
    dydt(3) = dr_dt;
    dydt(4) = R_precip;
    dydt(5) = Md_inflow;
    dydt(6) = Q*params.rho - E*params.rho;
end

function [value,isterminal,direction] = bot_closure_events(~,y,params)
    active = sum(y(3:6:end)>params.r_min);
    value = active-3;
    isterminal = 1;
    direction = -1;

    value = value(:)'; 
    isterminal = isterminal(:)'; 
    direction = direction(:)';
end
