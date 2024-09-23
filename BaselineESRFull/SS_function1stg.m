function g = SS_function1stg(~,x)
    
    Parameters % Load the parameters

    % CALCULATION ----------

    g = zeros(ns+1, 1); % Initialization of the function output
       
    % Reactions    R1 R2 R3 R4     Stoichiometric matrix           
    stoich_mat = [-1, -1, 0,  0;...  % C2H5OH
                          0,  0, -1, -3;... % H2O,
                          0,  1,  0,  0;...  % CH4
                          1,  1,  1,  5;...  %  H2,
                          0,  1, -1,  0;...  % CO
                          0,  0,  1,  2;...  % CO2
                          1,  0,  0, -1];  %  CH3CHO
    
    F = x(1:ns);  % Positions in vector x for molar flow rates 
    T = x(ns+1); % Position in vector x for temperature
    sumF = sum(F); % Total molar flow rate [mol/min]
    
    kreact = zeros(nr, 1);  % Kinetic parameters calculation
    for i = 1:nr
        kreact(i) = kinf(i) * exp(-Ea(i) * (1/(R*T) - 1/(R*Tref)));
    end
	kWGS = exp((4577.8/T) - 4.33);
   
   % Reaction rates determination
   r = zeros(nr, 1);
   Fratio_C2H5OH = F(1)/sumF; % Molar fractions
   Fratio_H2O = F(2)/sumF;
   Fratio_H2 = F(4)/sumF;
   Fratio_CO = F(5)/sumF;
   Fratio_CO2 = F(6)/sumF;
   Fratio_CH3CHO = F(7)/sumF;

   r(1) = kreact(1)*(P_bar*Fratio_C2H5OH); 
   r(2) = kreact(2)*(P_bar*Fratio_C2H5OH);
   r(3) = kreact(3)*((P_bar*Fratio_CO) *(P_bar*Fratio_H2O) - ...
                ((P_bar*Fratio_CO2)*(P_bar*Fratio_H2))/kWGS);
   r(4) = kreact(4)*(P_bar*Fratio_CH3CHO)*(P_bar*Fratio_H2O)^3;
            
   g(1:ns) = A*stoich_mat*r;  % Steady-state mass balance
   
   % Molar heat capacities at constant pressure % [J/(mol-K)]
   Cp_C2H5OH = 1.7690e1 + 1.4953e-1*T + 8.9481e-5*T^2 - 1.9738e-7*T^3 + 8.3175e-11*T^4;
   Cp_H2O = 3.4047e1 - 9.6506e-3*T + 3.2998e-5*T^2 - 2.0447e-8*T^3 + 4.3023e-12*T^4;
   Cp_CH4 = 3.8387e1 - 7.3663e-2*T + 2.9098e-4*T^2 - 2.6384e-7*T^3 + 8.0067e-11*T^4;
   Cp_H2 = 1.7639e1 + 6.7005e-2*T - 1.3148e-4*T^2 + 1.0588e-7*T^3 - 2.9180e-11*T^4;
   Cp_CO = 2.9006e1 + 2.4923e-3*T - 1.8644e-5*T^2 + 4.7989e-8*T^3 - 2.8726e-11*T^4;
   Cp_CO2 = 1.9022e1 + 7.9629e-2*T - 7.3706e-5*T^2 + 3.7457e-8*T^3 - 8.1330e-12*T^4;
   Cp_CH3CHO =2.4528e1 + 7.6013e-2*T + 1.3625e-4*T^2 - 1.9994e-7*T^3 + 7.5955e-11*T^4;
   Cp = [Cp_C2H5OH, Cp_H2O, Cp_CH4, Cp_H2, Cp_CO, Cp_CO2, Cp_CH3CHO]; 

   % Enthalpy of reactions at T
   deltaCp_rxn1 = (1/1)*Cp(4) + (1/1)*Cp(7) - (1)*Cp(1);
   deltaCp_rxn2 = (1/1)*Cp(4) + (1/1)*Cp(3) + (1/1)*Cp(5) - (1)*Cp(1);
   deltaCp_rxn3 = (1/1)*Cp(4) + (1/1)*Cp(6) - (1/1)*Cp(2) - (1)*Cp(5);
   deltaCp_rxn4 = (5/1)*Cp(4) + (2/1)*Cp(6) - (3/1)*Cp(2) - (1)*Cp(7);
   deltaCp = [deltaCp_rxn1, deltaCp_rxn2, deltaCp_rxn3, deltaCp_rxn4];
   deltaH = zeros(1, nr);
   for i = 1:nr
       deltaH(i) = deltaH_std(i) + deltaCp(i)*(T-298); % [J/mol] 
   end

   sumdFdz = sum(g(1:ns)); 
   sum3 = deltaH*r + sum(stoich_mat*r)*(-R*T);
   g(ns+1) = (U*a*(T_a-T)-((1/A)*R*T*sumdFdz)-sum3)/((1/A)*(Cp*F)); 
   % Steady-state energy balance

end