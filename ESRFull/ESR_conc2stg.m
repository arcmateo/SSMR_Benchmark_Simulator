function dxdt = ESR_conc2stg(~,x,x_boundary,p,F_in)

    p = Parameters(); % Load the parameters 
    np = p.np; % Number of points (spatial discretization)
    R = p.R; % [J/(mol-K)] Gas universal constant
    Patm =p.Patm; % [Pa] Atmosferic pressure   
    P_in = p.P_in; % [Pa] Inlet pressure
    T_a = p.T_a; % [K] Furnace temperature 
    U = p.U; % [J /(m2-min-K)] Heat transfer coefficient   
    A = p.A; % [m2] Reactor cross-sectional area      
    a = p.a; % [m2/m3] % area per reactor volume for heat transfer    
    ns = p.ns; % Number of species  
    deltaz1= p.deltaz1; % delta_z 1st stage   
    deltaz2 = p.deltaz2; % delta_z 2nd stage   
    pe0 = p.pe0; % [mol / (min-Pa^1/2-m)] Pre-exponential factor   
    Eam = p.Eam; % [J/mol] "Activation energy" of membrane   
    deltam = p.deltam; % [m] membrane thickness   
    Dm = p.Dm; % [m] Membrane diameter

   
    % CALCULATION

    C_in = x_boundary(1:ns); % [mol/m3] Vector of inlet Concentrations
    T_in = x_boundary(ns+1); % [K] Inlet temperature

    v_in = ((R*T_in)/(A*P_in))*sum(F_in); % [m/min] Inlet velocity

    v_k_0 = 0; % velocity must be initialized to something to record v_k_1 for k > 1

    dxdt = zeros((ns+1)*np, 1);

    for k = 1:np 
       C_k = x(k:np:ns*np); % Vector of conc. of each species at point k
       index3 = ns*np + k; % Index for the temperature at point k
       T_k = x(index3); % Temperature at point k

       pe_k = pe0 * exp(-Eam/(R*T_k)); % [mol/m min Pa^(1/2)]
       P_H2_k = C_k(4)*R*T_k; % [Pa]
       
       F_H2_perm_k = pe_k*((pi*Dm*deltaz2)/deltam)*...
           (sqrt(P_H2_k)-sqrt(Patm)); % [mol/min]

       if F_H2_perm_k < 0
          F_H2_perm_k = 0;
       end
       
       F_H2_perm_k_vol = (F_H2_perm_k)/(A*deltaz2);

       if k == 1
           Ftot_k = sum(F_in) - F_H2_perm_k;
       else
           Ftot_k = Ftot_k - F_H2_perm_k;
       end
   
       v_k_1 = v_k_0; % velocity at previous point must be recorded
       v_k_0 = ((R*T_k)/(P_in*A))*Ftot_k; % [m/min] v = (RT/PA)*Flow
       
       % Computing dCdt
       for j = 1:ns % for species j
           index = (j-1)*np + k;
           if k == 1
               dvC_dz_j_k = (v_k_0*x(index) - v_in*C_in(j)) / deltaz2;
           else
               dvC_dz_j_k = (v_k_0*x(index) - v_k_1*x(index - 1)) / deltaz2;
           end

           if j == 4
               dxdt(index) = -dvC_dz_j_k - F_H2_perm_k_vol;
           else
               dxdt(index) = -dvC_dz_j_k;
           end

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

       % Molar internal energy of Hydrogen
       U_H2 = Cv(4)*(T_k - 0);
  
       % Computing dTdt
   
       if k == 1
           dTdz_k = (x(index3) - T_in) / deltaz1;
       else
           dTdz_k = (x(index3) - x(index3 - 1)) / deltaz1;
       end  
       dxdt(index3) = (U*a*(T_a-T_k) - (Cp*C_k)*v_k_0*dTdz_k...
           - F_H2_perm_k_vol*U_H2) / (Cv*C_k);    
   end

end
