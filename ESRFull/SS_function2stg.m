function g = SS_function2stg(~,x,p)
 
    R = p.R; % [J/(mol-K)] Gas universal constant
    Patm =p.Patm; % [Pa] Atmosferic pressure
    P_in = p.P_in; % [Pa] Inlet pressure    
    T_a = p.T_a; % [K] Furnace temperature   
    U = p.U; % [J /(m2-min-K)] Heat transfer coefficient
    A = p.A; % [m2] Reactor cross-sectional area
    a = p.a; % [m2/m3] % area per reactor volume for heat transfer
    ns = p.ns; % Number of species
    pe0 = p.pe0; % [mol / (min-Pa^1/2-m)] Pre-exponential factor
    Eam = p.Eam; % [J/mol] "Activation energy" of membrane
    deltam = p.deltam; % [m] membrane thickness
    Dm = p.Dm; % [m] Membrane diameter


    % CALCULATION ----------

    g = zeros(ns+1, 1); % Initialization of the function output
    
    F = x(1:ns); % [mol/min] Positions for molar flows
    T = x(ns+1); % [K] Position for temperature
    sumF = sum(F); % [mol/min] Total molar flow
    
    H2 = 4; % Hydrogen index according to the order established here
    pe = pe0 * exp(-Eam / (R*T)); % [[mol/(min-Pa^1/2-m)]] Permeability
    F_H2_perm_deltaz = (pe/deltam)*(pi*Dm)*(sqrt(P_in*F(H2)/sumF) - sqrt(Patm));
    % [mol/min-m] Flow of hydrogen modeled by the Sievert's law

   % Molar heat capacities at constant pressure % [J/(mol-K)]
   Cp_C2H5OH = 1.7690e1 + 1.4953e-1*T + 8.9481e-5*T^2 - 1.9738e-7*T^3 + 8.3175e-11*T^4;
   Cp_H2O = 3.4047e1 - 9.6506e-3*T + 3.2998e-5*T^2 - 2.0447e-8*T^3 + 4.3023e-12*T^4;
   Cp_CH4 = 3.8387e1 - 7.3663e-2*T + 2.9098e-4*T^2 - 2.6384e-7*T^3 + 8.0067e-11*T^4;
   Cp_H2 = 1.7639e1 + 6.7005e-2*T - 1.3148e-4*T^2 + 1.0588e-7*T^3 - 2.9180e-11*T^4;
   Cp_CO = 2.9006e1 + 2.4923e-3*T - 1.8644e-5*T^2 + 4.7989e-8*T^3 - 2.8726e-11*T^4;
   Cp_CO2 = 1.9022e1 + 7.9629e-2*T - 7.3706e-5*T^2 + 3.7457e-8*T^3 - 8.1330e-12*T^4;
   Cp_CH3CHO =2.4528e1 + 7.6013e-2*T + 1.3625e-4*T^2 - 1.9994e-7*T^3 + 7.5955e-11*T^4;
   Cp = [Cp_C2H5OH, Cp_H2O, Cp_CH4, Cp_H2, Cp_CO, Cp_CO2, Cp_CH3CHO];
                
    if F_H2_perm_deltaz < 0
        F_H2_perm_deltaz = 0; % It is not possible to "suck" H2 from the output
    end
    
    g(H2) = (-F_H2_perm_deltaz); % Mass Balance SS
    
    g(ns+1) = (U*a*(T_a - T)) / ((1/A)*(Cp* F)); % Energy Balance SS
                
end
