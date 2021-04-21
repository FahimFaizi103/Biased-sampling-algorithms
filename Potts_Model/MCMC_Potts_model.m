function[Mmean, Emean, spin] = MCMC_Potts_model(numSweeps, spin, kT, J, q_states, K, X_equilibrium)

Mmean = zeros();
Emean = zeros();

for SweepIndex = 1:1:numSweeps
    
    for i = 1:1:numel(spin)
        
        sample = i;
        
        % LinearIndex picks a spin by acing as an index into the Matrix: spin. We
        % can convert it into [row, col] of our matrix spin by the following:
        
        [row,col] = ind2sub(size(spin),sample);
        
        % The subscripts of the neighbours of spin(row,col) are the following.
        % The "above" and "below" identify the rows that the nearest neighbours of
        %spin(row,col) lies, likewise the "left" and "right" identify the coloums in
        %which the nearest neighbours lie.
        
        above = mod(row - 1 - 1, size(spin, 1)) + 1;
        below = mod(row + 1 - 1, size(spin, 1)) + 1;
        left  = mod(col - 1 - 1, size(spin, 2)) + 1;
        right = mod(col + 1 - 1, size(spin, 2)) + 1;
        
        % We now specify the state of the nearest neighbours.
        neighbors = [       spin(above, col);
            spin(row, left);                 spin(row, right);
            spin(below, col)];
        
        
        [Boltzmann_weights] = Boltzmann_weights_Potts(q_states, neighbors, kT, J, spin, row, col, K, X_equilibrium);
            
        P =  Boltzmann_weights./(sum(Boltzmann_weights));
        
        new_state = q_states(find(rand<cumsum(P),1,'first'));
        
        spin(row,col) = new_state;
        
        
    end
    
    Mmean(SweepIndex) = mean(spin(:)); % compute magnetisation per spin.
    [Emean(SweepIndex)] = Potts_energy_compute(spin,J);
    
    
    
end


end

