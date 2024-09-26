function f = ESR_flows2stg(~,x,x_boundary)
    
    p = Parameters(); % Load the parameters 
    np = p.np; % Number of points (spatial discretization)
    R = p.R; % [J/(mol-K)] Gas universal constant
    Patm =p.Patm; % [Pa] Atmosferic pressure
    P0 = p.P0; % [Pa] Inlet pressure  
    T_a = p.T_a; % [K] Furnace temperature
    U = p.U; % [J /(m2-min-K)] Heat transfer coefficient
    A = p.A; % [m2] Reactor cross-sectional area  
    a = p.a; % [m2/m3] % area per reactor volume for heat transfer
    ns = p.ns; % Number of species
    deltaz2 = p.deltaz2; % delta_z 2nd stage
    pe0 = p.pe0; % [mol / (min-Pa^1/2-m)] Pre-exponential factor
    Eam = p.Eam; % [J/mol] "Activation energy" of membrane
    deltam = p.deltam; % [m] membrane thickness
    Dm = p.Dm; % [m] Membrane diameter


    % CALCULATION
    
    f = zeros((ns+1)*np, 1);
    
    F_in = x_boundary(1:ns);
    T_in = x_boundary(ns+1);
    
    for k = 1 : np
        F_k = x(k:np:ns*np);
        sumF_k = sum(F_k);
        index3 = ns*np + k;
        T_k = x(index3);
        
        pe_k = pe0 * exp(-Eam / (R*T_k));
    	F_H2_k = x((4-1)*np + k);
        %if (F_H2_k < 0 && abs(F_H2_k) < 1e-10)
        %    F_H2_k = 0; % Having problems with -6.3e-13 in F_H2
        %end
    	F_H2_perm_deltaz_k = (pe_k / deltam) * (pi * Dm) * ...
                        (sqrt(P0 * F_H2_k/sumF_k) - sqrt(Patm));
        %if F_H2_perm_deltaz_k < 0 % TEMP TEMP TEMP
        %    F_H2_perm_deltaz_k = 0;
        %end
        
        sumdFdz_k = 0;
        for j = 1 : ns
            index = (j-1)*np + k;
            if k == 1
                dFdz_j_k = (x(index) - F_in(j)) / deltaz2;
            else
                dFdz_j_k = (x(index) - x(index - 1)) / deltaz2;
            end
            
            if j == 4
                f(index) = R * T_k * sumF_k * ...
                    	1/A * (- F_H2_perm_deltaz_k - dFdz_j_k);
            else
                f(index) = R * T_k * sumF_k * (- 1/A * dFdz_j_k);
            end
            
            sumdFdz_k = sumdFdz_k + dFdz_j_k;
        end

        if k == 1
            dTdz_k = (x(index3) - T_in) / deltaz2;
        else
            dTdz_k = (x(index3) - x(index3 - 1)) / deltaz2;
        end

        % Molar heat capacities at constant pressure % [J/(mol-K)]
        Cp_C2H5OH = 1.7690e1 + 1.4953e-1*T_k + 8.9481e-5*T_k^2 - 1.9738e-7*T_k^3 + 8.3175e-11*T_k^4;
        Cp_H2O = 3.4047e1 - 9.6506e-3*T_k + 3.2998e-5*T_k^2 - 2.0447e-8*T_k^3 + 4.3023e-12*T_k^4;
        Cp_CH4 = 3.8387e1 - 7.3663e-2*T_k + 2.9098e-4*T_k^2 - 2.6384e-7*T_k^3 + 8.0067e-11*T_k^4;
        Cp_H2 = 1.7639e1 + 6.7005e-2*T_k - 1.3148e-4*T_k^2 + 1.0588e-7*T_k^3 - 2.9180e-11*T_k^4;
        Cp_CO = 2.9006e1 + 2.4923e-3*T_k - 1.8644e-5*T_k^2 + 4.7989e-8*T_k^3 - 2.8726e-11*T_k^4;
        Cp_CO2 = 1.9022e1 + 7.9629e-2*T_k - 7.3706e-5*T_k^2 + 3.7457e-8*T_k^3 - 8.1330e-12*T_k^4;
        Cp_CH3CHO =2.4528e1 + 7.6013e-2*T_k + 1.3625e-4*T_k^2 - 1.9994e-7*T_k^3 + 7.5955e-11*T_k^4;
        Cp = [Cp_C2H5OH, Cp_H2O, Cp_CH4, Cp_H2, Cp_CO, Cp_CO2, Cp_CH3CHO]; 
        
        f(index3) = U * a * (T_a - T_k) - 1/A * R * T_k * sumdFdz_k - ...
                        1/A * R * T_k * F_H2_perm_deltaz_k - ...
                        1/A * (Cp * F_k) * dTdz_k;
    end
    
end
