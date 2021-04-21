
clear all
clc


J = 1;
kT = 2.0;
numSim = 1;
numSpinsPerDim = 12;
probSpinUp = 0.5;
numSweeps = 10^3;

spin = ones(numSpinsPerDim, numSpinsPerDim).*-1;

%% Biasing parameters

num_windows = 10;
lower_boundary = -1;
upper_boundary =  1;
X_equilibrium = linspace(lower_boundary, upper_boundary, num_windows);
X_equilibrium_record = zeros();
K = ones(1,num_windows).*200; % Hook's constant.
%% Biased sampling
Mmean = zeros(num_windows, numSweeps);
Emean = zeros(num_windows, numSweeps);

for i = 1:1:num_windows
    
   
    [Mmean(i,:), Emean(i,:), spin] = MCMC_ising_model(kT, numSweeps, J, spin, K(i), X_equilibrium(i));
    
    X_equilibrium_record(i) = X_equilibrium(i);
    
    
end

