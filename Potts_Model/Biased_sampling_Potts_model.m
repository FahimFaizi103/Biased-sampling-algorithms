
clear all
clc

J = 1;
kT = 0.8;
numSim = 1;
numSpinsPerDim = 5;
numSweeps = 10^3;

min_q = 1;
max_q = 4;

spin = ones(numSpinsPerDim,numSpinsPerDim);
q_states = (min_q:1:max_q);

%% Biasing parameters

num_windows = 100;
lower_boundary =  1;
upper_boundary =  4;
X_equilibrium = linspace(lower_boundary, upper_boundary, num_windows);
X_equilibrium_record = zeros();
K = ones(1,num_windows).*100; % Hook's constant.
%%
Mmean = zeros(num_windows, numSweeps);
Emean = zeros(num_windows, numSweeps);

for i = 1:1:num_windows
    
   
    [Mmean(i,:), Emean(i,:), spin] = MCMC_Potts_model(numSweeps, spin, kT, J, q_states, K(i), X_equilibrium(i));
    
    X_equilibrium_record(i) = X_equilibrium(i);
    
    
end

