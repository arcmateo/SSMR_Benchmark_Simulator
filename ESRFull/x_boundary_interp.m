function x_interp = x_boundary_interp(x_boundary, t_boundary, t)

    %%% NOTE: Each row is a full state vector for a certain time instant!!
    num_states = size(x_boundary,2);
    x_interp = zeros(1,num_states);
    for state = 1 : num_states
        state_boundary = x_boundary(:,state);
        x_interp(state) = interp1(t_boundary, state_boundary, t);
    end
end