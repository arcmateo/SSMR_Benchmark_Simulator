function p = Parameters()

% PARAMETERS ----------

% Order of reactions: [dehydrogenation, decomposition, WGSR, reforming]
% Order of species: [C2H5OH, H2O, CH4, H2, CO, CO2, CH3CHO]

p.np = 50; % Number of points (spatial discretization)

p.R = 8.31432; % [J/(mol-K)] Gas universal constant

p.Patm = 101325; % [Pa] Atmosferic pressure

p.P_in = 4*p.Patm; % [Pa] Inlet pressure

p.P_bar = p.P_in/p.Patm; % [bar] Constant pressure assumed along the reactor
    
p.T_in = 873.15; % [K] Inlet temperature
    
p.Tref = 773.15; % [K] Reference temperature
    
p.T_a = p.T_in; % [K] Furnace temperature
    
p.Ea = [7.0e4, 1.30e5, 7.0e4, 9.8e4]; % [J/mol] Activation energy
    
p.kinf = [2.1e4, 2.0e3, 1.9e4, 2.0e5]; % [mol/(m3-min-bar)]
% Pre-exponential factor
    
p.deltaH_std = [64600, 49875, -41166, 109136]; % [J/mol] 
% Standard enthalpy of reactions

p.U = (25)*60; % [J /(m2-min-K)] Heat transfer coefficient

p.d = 22e-3; % [m] Reactor diameter

p.A = pi*((p.d^2)/4); % [m2] Reactor cross-sectional area
    
p.a = 4/p.d; % [m2/m3] % area per reactor volume for heat transfer

p.ns = 7; % Number of species

p.nr = 4; % Number of reactions

p.L = 0.23; % [m] Reactor length

p.L2 = 0.076; % [m] 2nd stage length

p.L1 = p.L - p.L2; % [m] 1st stage length

p.deltaz1 = p.L1 / p.np; % delta_z 1st stage

p.deltaz2 = p.L2 / p.np; % delta_z 2nd stage

p.pe0 = (2.25e-8)*60; % [mol / (min-Pa^1/2-m)] Pre-exponential factor

p.Eam = 8.8e3; % [J/mol] "Activation energy" of membrane

p.deltam = 3.0e-5; % [m] membrane thickness

p.Dm = ((1/8)*0.0254); % [m] Membrane diameter

end