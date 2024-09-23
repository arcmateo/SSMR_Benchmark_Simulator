% PARAMETERS ----------

% Order of reactions: [dehydrogenation, decomposition, WGSR, reforming]
% Order of species: [C2H5OH, H2O, CH4, H2, CO, CO2, CH3CHO]

np = 600; % Number of points (spatial discretization)

R = 8.31432; % [J/(mol-K)] Gas universal constant

Patm = 101325; % [Pa] Atmosferic pressure

P0 = 4*Patm; % [Pa] Inlet pressure

P_bar = P0/Patm; % [bar] Constant pressure assumed along the reactor
    
T0 = 873.15; % [K] Inlet temperature
    
Tref = 773.15; % [K] Reference temperature
    
T_a = T0; % [K] Furnace temperature
    
Ea = [7.0e4, 1.30e5, 7.0e4, 9.8e4]; % [J/mol] Activation energy
    
kinf = [2.1e4, 2.0e3, 1.9e4, 2.0e5]; % [mol/(m3-min-bar)]
% Pre-exponential factor
    
deltaH_std = [64600, 49875, -41166, 109136]; % [J/mol] 
% Standard enthalpy of reactions

U = (25)*60; % [J /(m2-min-K)] Heat transfer coefficient

d = 22e-3; % [m] Reactor diameter

A = pi*((d^2)/4); % [m2] Reactor cross-sectional area
    
a = 4/d; % [m2/m3] % area per reactor volume for heat transfer

T0 = 873.15; % [K] Inlet temperature

ns = 7; % Number of species

nr = 4; % Number of reactions

L = 0.23; % [m] Reactor length

L2 = 0.076; % [m] 2nd stage length

L1 = L - L2; % [m] 1st stage length

deltaz1 = L1 / np; % delta_z 1st stage

deltaz2 = L2 / np; % delta_z 2nd stage

pe0 = (2.25e-8)*60; % [mol / (min-Pa^1/2-m)] Pre-exponential factor

Eam = 8.8e3; % [J/mol] "Activation energy" of membrane

deltam = 3.0e-5; % [m] membrane thickness

Dm = ((1/8)*0.0254); % [m] Membrane diameter