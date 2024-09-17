function xc = Flows_to_Conc(xf)
% Transforms flows to concentrations as a ROW VECTOR
%  xf: vector of flows 
%  xc: vector of concentrations 

% PARAMETERS ----------

np = 600; % Number of points (spatial discretization)

R = 8.31432; % [J/(mol-K)] Gas universal constant

Patm = 101325; % [Pa] Atmosferic pressure

P0 = 4*Patm; % [Pa] Inlet pressure

ns = 7; % Number of species

xc = zeros(1,np*(ns+1)); % initialize the vector of concentrations and temperature

for k = 1:np % at every discretization point k
    F_k = xf(k:np:ns*np); % F_k is a vector of flows of each species at point k
    sumF_k = sum(F_k); % Total molar flow at point k
    index3 = ns*np + k; % Index for the temperature at point k
    T_k = xf(index3); % Temperature at point k
    xc(k:np:ns*np) = ((F_k/sumF_k)/(R*T_k))*P0;  % [mol/m3]
    xc(index3) = T_k;
end

end

 


