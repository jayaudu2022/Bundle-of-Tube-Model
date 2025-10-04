% ---Single capillary tube ---
clear; clc;

%% 1. Parameters
case_name = 'radius_75um';

pars.T      = 323.15;       % Temperature (K)
pars.Pg     = 101325;       % Gas pressure (Pa)
pars.RH     = 0.0;          % Relative humidity 
pars.rho    = 1000;         % Brine density (kg/m3)
pars.sigma  = 0.0718;       % Surface tension (N/m)
pars.mu     = 1.0142e-3;    % Viscosity (Pa·s)
pars.costh  = cosd(30);     % Contact angle 
pars.g      = 9.81;         % Gravity (m/s2)
pars.rho_s  = 2165;         % Salt density (kg/m3)
pars.Vm     = 1.802e-5;     % Molar volume (m3/mol)
pars.Rg     = 8.314;        % Universal Gas constant (J/mol/K)
pars.Psat   = 12.33e3;      % Saturation vapor pressure (Pa)
pars.Mw     = 0.018;        % Molar mass of water (kg/mol)
pars.Dv     = 2.5e-4;       % Diffusivity of vapor (m2/s)
pars.L      = 3.0;          % Tube length (m)

% Salt properties (molality-based)
pars.m0     = 6.0;          % Inflow molality (mol/kg)
pars.m_sat  = 6.14;         % Saturation molality (mol/kg)
pars.M_NaCl = 0.05844;      % Molar mass of NaCl (kg/mol)

% Geometry
pars.h0     = 5e-3;         % Initial meniscus height (m)
pars.r0     = 75e-6;        % Initial radius (m)
pars.r_min  = 5e-6;         % Minimum radius (m)

% Precipitation parameters
pars.tau_inst    = 1e-9;
pars.l_rim       = 1e-3;
pars.z_min       = 1e-6;
pars.use_gravity = 1;

% Output settings
pars.output_interval = 50000;
pars.plot_interval   = 100;

%% 2. Initial Conditions
V0       = pi * pars.r0^2 * pars.h0;
M_water0 = pars.rho * V0;
Md0      = pars.m0 * M_water0 * pars.M_NaCl;

y0 = [pars.h0; Md0; pars.r0; 0; 0; M_water0];
tspan = [0, 1e8];
options = odeset('RelTol', 1e-5, 'MaxStep', 1, ...
    'Events', @(t,y) radius_drop_event(t,y,pars));

%% 3. Run Simulation
fprintf('=== SALT PRECIPITATION SIMULATION ===\n');
fprintf('Initial molality: %.2f (%.1f%% of sat)\n', pars.m0, pars.m0/pars.m_sat*100);

[t, y] = ode15s(@(t,y) capillaryTransportODE(t,y,pars), tspan, y0, options);

% Extract results
h_all      = y(:,1);
Md_all     = y(:,2);
r_all      = y(:,3);
Mprecip_all = y(:,4);
M_inflow_all = y(:,5);
M_water_all  = y(:,6);

molality_all = (Md_all ./ pars.M_NaCl) ./ max(M_water_all, 1e-20);

% Mass balance check
M_system_total = Md_all(end) + Mprecip_all(end);
M_inflow_total = M_inflow_all(end) + Md0;
M_balance_error = abs(M_system_total - M_inflow_total) ...
                  / max(M_inflow_total,1e-20) * 100;

fprintf("Final mass balance error = %.2f%%\n", M_balance_error);

save(['results_' case_name '.mat'], ...
    't','y','h_all','r_all','Md_all','Mprecip_all', ...
    'M_inflow_all','M_water_all','molality_all','pars','case_name');

%% ----------------------- FUNCTIONS -----------------------

function dydt = capillaryTransportODE(t, y, pars)
    % Core transport model: imbibition, evaporation, precipitation
    
    % Unpack state
    h = max(y(1),0);
    Md = max(y(2),0);
    r = max(y(3),pars.r_min);
    M_water = y(6);

    dydt = zeros(6,1);

    % Geometry
    h_eff = min(h, 0.9*pars.L);
    A = pi * r^2;
    V = A * h_eff;

    % Capillary pressure and resistance
    Pc = 2*pars.sigma*pars.costh / r;
    pgrav = pars.use_gravity * pars.rho * pars.g * h_eff;
    Rh = 8*pars.mu*h_eff / (pi*r^4);

    % Evaporation (Stefan diffusion + Kelvin correction)
    kelv = min((2*pars.sigma*pars.Vm*pars.costh)/(r*pars.Rg*pars.T), 50);
    y1 = min(pars.RH*pars.Psat/pars.Pg, 0.999999);
    y2 = min(pars.Psat*exp(kelv)/pars.Pg, 0.999999);
    lnterm = max(log((1-y1)/(1-y2)), 0);
    L_eff = max(pars.L - h_eff, 1e-9);
    Jm = (pars.Pg*pars.Dv*pars.Mw)/(pars.Rg*pars.T) * (lnterm/L_eff);
    E = Jm*A/pars.rho;   % (m3/s)

    % Flow rate & height ODE
    Q = (Pc - pgrav)/max(Rh,1e-30);
    dydt(1) = (Q - E)/A;

    % Salt inflow
    Md_inflow = pars.m0 * Q * pars.rho * pars.M_NaCl;

    % Precipitation
    mol_current = (Md/pars.M_NaCl)/max(M_water,1e-20);
    if mol_current > pars.m_sat
        supersat = mol_current - pars.m_sat;
        k_precip = 1e-4;
        if mol_current > 1.5*pars.m_sat, k_precip = k_precip*100;
        elseif mol_current > 1.2*pars.m_sat, k_precip = k_precip*10;
        end
        R_precip = min(k_precip*supersat*M_water*pars.M_NaCl, 0.9*Md);
    else
        R_precip = 0;
    end

    % Radius change from precipitation
    if R_precip>0
        z_dep = min(max(h_eff,0), max(pars.z_min, pars.l_rim));
        dr_dt = -R_precip/(pars.rho_s*2*pi*r*z_dep);
    else
        dr_dt = 0;
    end

    % Assemble ODEs
    dydt(2) = Md_inflow - R_precip;          % Salt mass in solution
    dydt(3) = dr_dt;                         % Radius change
    dydt(4) = R_precip;                      % Precipitated salt
    dydt(5) = Md_inflow;                     % Cumulative inflow
    dydt(6) = Q*pars.rho - E*pars.rho;       % Water balance
end

function [value,isterminal,direction] = radius_drop_event(~,y,pars)
    r = y(3);
    value = r - pars.r_min;
    isterminal = 1;
    direction = -1;
end
