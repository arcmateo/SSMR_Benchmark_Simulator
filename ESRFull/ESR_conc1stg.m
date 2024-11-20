function dxdt = ESR_conc1stg(t,x,u,p)

    global F_out

    np = p.np; % Number of points (spatial discretization)
    R = p.R; % [J/(mol-K)] Gas universal constant    
    Patm =p.Patm; % [Pa] Atmosferic pressure    
    P_in = p.P_in; % [Pa] Inlet pressure       
    Tref = p.Tref; % [K] Reference temperature        
    T_a = p.T_a; % [K] Furnace temperature        
    Ea = p.Ea; % [J/mol] Activation energy       
    kinf = p.kinf; % [mol/(m3-min-bar)]
    % Pre-exponential factor       
    deltaH_std = p.deltaH_std; % [J/mol] 
    % Standard enthalpy of reactions   
    U = p.U; % [J /(m2-min-K)] Heat transfer coefficient    
    A = p.A; % [m2] Reactor cross-sectional area   
    a = p.a; % [m2/m3] % area per reactor volume for heat transfer
    T_in = p.T_in; % [K] Inlet temperature
    ns = p.ns; % Number of species
    nr = p.nr; % Number of reactions
    deltaz1= p.deltaz1; % delta_z 1st stage
    
    % Input change (change to simulate a manipulation)
    if t > 0.5
       u(1) = 0.0021;
    end

    F_in = [u(1), u(2), zeros(1,ns-2)]; % [mol/min] Vector of inlet molar flow rates
    C_in = ((F_in/sum(F_in))/(R*T_in))*P_in; % [mol/m3] Vector of inlet Concentrations
    v_in = ((R*T_in)/(A*P_in))*sum(F_in); % [m/min] Inlet velocity

    % Reactions    R1 R2 R3 R4     Stoichiometric matrix           
    stoich_mat = [-1, -1, 0,  0;...  % C2H5OH
                          0,  0, -1, -3;... % H2O,
                          0,  1,  0,  0;...  % CH4
                          1,  1,  1,  5;...  %  H2,
                          0,  1, -1,  0;...  % CO
                          0,  0,  1,  2;...  % CO2
                          1,  0,  0, -1];  %  CH3CHO

    v_k_0 = 0; % velocity must be initialized to something to record v_k_1 for k > 1

    dxdt = zeros((ns+1)*np, 1);

    for k = 1:np 
       C_k = x(k:np:ns*np); % Vector of conc. of each species at point k
       index3 = ns*np + k; % Index for the temperature at point k
       T_k = x(index3); % Temperature at point k
       
       kreact_k = zeros(nr, 1); % vector of kinetics parameters
       for i = 1:nr
           kreact_k(i) = kinf(i)*exp(-Ea(i)*(1/(R*T_k)-1/(R*Tref)));
       end
       kWGS_k = exp((4577.8/T_k) - 4.33);
       
       r_k = zeros(nr, 1); % Reaction rates at point k

       % Partial pressure of each especies in [bar]
       P_C2H5OH = (C_k(1)*R*T_k)/Patm;
       P_H2O = (C_k(2)*R*T_k)/Patm;
       % Partial pressure of methane is not required
       P_H2 = (C_k(4)*R*T_k)/Patm;
       P_CO = (C_k(5)*R*T_k)/Patm;
       P_CO2 = (C_k(6)*R*T_k)/Patm;
       P_CH3CHO = (C_k(7)*R*T_k)/Patm;
   
       r_k(1) = kreact_k(1)*(P_C2H5OH); 
       r_k(2) = kreact_k(2)*(P_C2H5OH); 
       r_k(3) = kreact_k(3)*(P_CO*P_H2O - ((P_CO2*P_H2)/kWGS_k));
       r_k(4) = kreact_k(4)*P_CH3CHO*P_H2O^3;
       
       molarity_change = r_k(1) + 2*r_k(2) + 3*r_k(4);
   
       if k == 1
           Ftot_k = sum(F_in) + A*molarity_change*deltaz1;
       else
           Ftot_k = Ftot_k + A*molarity_change*deltaz1;
       end
        
       v_k_1 = v_k_0; % velocity at previous point must be recorded
       v_k_0 = ((R*T_k)/(P_in*A))*Ftot_k; % [m/min] v = (RT/PA)*Flow
       
       % Computing dCdt
       for j = 1:ns % for species j
           index = (j-1)*np + k;
           if k == 1
               dvC_dz_j_k = (v_k_0*x(index) - v_in*C_in(j)) / deltaz1;
           else
               dvC_dz_j_k = (v_k_0*x(index) - v_k_1*x(index - 1)) / deltaz1;
           end
   
           rj_k = stoich_mat(j,:)*r_k;
   
           dxdt(index) = rj_k - dvC_dz_j_k;
       end
      
       % Update of properties that change with temperature
      
       % Molar heat capacities at constant pressure % [J/(mol-K)]
       Cp_C2H5OH = 1.7690e1 + 1.4953e-1*T_k + 8.9481e-5*T_k^2 - 1.9738e-7*T_k^3 + 8.3175e-11*T_k^4;
       Cp_H2O = 3.4047e1 - 9.6506e-3*T_k + 3.2998e-5*T_k^2 - 2.0447e-8*T_k^3 + 4.3023e-12*T_k^4;
       Cp_CH4 = 3.8387e1 - 7.3663e-2*T_k + 2.9098e-4*T_k^2 - 2.6384e-7*T_k^3 + 8.0067e-11*T_k^4;
       Cp_H2 = 1.7639e1 + 6.7005e-2*T_k - 1.3148e-4*T_k^2 + 1.0588e-7*T_k^3 - 2.9180e-11*T_k^4;
       Cp_CO = 2.9006e1 + 2.4923e-3*T_k - 1.8644e-5*T_k^2 + 4.7989e-8*T_k^3 - 2.8726e-11*T_k^4;
       Cp_CO2 = 1.9022e1 + 7.9629e-2*T_k - 7.3706e-5*T_k^2 + 3.7457e-8*T_k^3 - 8.1330e-12*T_k^4;
       Cp_CH3CHO =2.4528e1 + 7.6013e-2*T_k + 1.3625e-4*T_k^2 - 1.9994e-7*T_k^3 + 7.5955e-11*T_k^4;
       Cp = [Cp_C2H5OH, Cp_H2O, Cp_CH4, Cp_H2, Cp_CO, Cp_CO2, Cp_CH3CHO]; 
       Cv = Cp - R;
       
       % Enthalpy of reactions 
       deltaCp_rxn1 = (1/1)*Cp(4) + (1/1)*Cp(7) - (1)*Cp(1);
       deltaCp_rxn2 = (1/1)*Cp(4) + (1/1)*Cp(3) + (1/1)*Cp(5) - (1)*Cp(1);
       deltaCp_rxn3 = (1/1)*Cp(4) + (1/1)*Cp(6) - (1/1)*Cp(2) - (1)*Cp(5);
       deltaCp_rxn4 = (5/1)*Cp(4) + (2/1)*Cp(6) - (3/1)*Cp(2) - (1)*Cp(7);
       deltaCp = [deltaCp_rxn1, deltaCp_rxn2, deltaCp_rxn3, deltaCp_rxn4];
       deltaH = zeros(1, nr);
       for i = 1:nr
          deltaH(i) = deltaH_std(i) + deltaCp(i)*(T_k-298); % [J/mol] 
       % deltaH = deltaH_std + deltaCp*(T-Tref)
       end
 
       % Computing dTdt
       Hr_k = -deltaH_std*r_k; % Total enthalpy of all reactions at point k
   
       if k == 1
           dTdz_k = (x(index3) - T_in) / deltaz1;
       else
           dTdz_k = (x(index3) - x(index3 - 1)) / deltaz1;
       end  
       dxdt(index3) = (U*a*(T_a-T_k) + Hr_k - (Cp*C_k) * v_k_0*dTdz_k) / (Cv*C_k);    
   end
   F_out = Ftot_k;
end
