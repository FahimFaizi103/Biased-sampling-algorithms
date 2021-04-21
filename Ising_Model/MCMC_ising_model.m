function[Mmean, Emean, spin] = MCMC_ising_model(kT, numSweeps, J, spin, K, X_equilibrium)
%--------------------------------------------------------------------------
% Sequential_Metropolis
%--------------------------------------------------------------------------
Mmean = zeros();
Emean = zeros();
%% Ising model MCMC simulations

for SweepIndex = 1:1:numSweeps
    
    for i = 1:1:numel(spin)
        
        
        [dE, row, col] = sample_spins_sequentially_biased(spin,J,i,K,X_equilibrium);
        
        [spin] = spin_flip_Met(dE, kT, spin, row, col);
        
    end
    
     Mmean(SweepIndex) = mean(spin(:));
    [Emean(SweepIndex)] = energy_compute(spin,J); % compute energy per spin.
    
    
end
