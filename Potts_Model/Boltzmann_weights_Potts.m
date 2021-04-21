function[Boltzmann_weights] = Boltzmann_weights_Potts(q_states, neighbors, kT, J, spin, row, col, K, X_equilibrium)

Mmean_ROL = mean(spin(:)) - spin(row,col)./(numel(spin)); % Mmean of the rest of the lattice.
Boltzmann_weights = zeros();

for m = 1:1:numel(q_states)
    
    state = q_states(m); 
    interaction_after_flip = neighbors;  
    
    interaction_after_flip(interaction_after_flip ~= state) = 0; % unlike neighbours contribute zeros to Hamiltonian
    interaction_after_flip(interaction_after_flip == state) = 1; % Like neighbours contribute -J to Hamiltonian.

    E_bath = -J.*(sum(interaction_after_flip)) + K.*(Mmean_ROL + state./numel(spin) - X_equilibrium).^2;
    Boltzmann_weights(m) = exp(-E_bath./kT);
     
    
end

end