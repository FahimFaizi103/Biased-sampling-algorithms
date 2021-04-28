%{
 Biased_sampling_Ising_model with addition of Dancing Umbrellas for the
 positions
%}
addpath(genpath("./"))
%%
% clear all
clc

J = 1;
kT = 2.0;
numSpinsPerDim = 12;
probSpinUp = 0.5;
numSweeps = 10^3;
% numSweeps = 2*10^2;
spin = ones(numSpinsPerDim, numSpinsPerDim).*-1;

%% Biasing parameters

num_targets = 6; % replaces previous num_windows
lower_boundary = -1;
upper_boundary =  1;
K_const = 50;

numUpdates =1;

%% Instantiation
num_windows = num_targets*(1+numUpdates);
numSweeps_per_split = floor(numSweeps./(1+numUpdates));

% X_equilibrium = linspace(lower_boundary, upper_boundary, num_windows);
Target_centers = linspace(lower_boundary+0.1, upper_boundary-0.1, num_targets);

K = ones(1,num_windows).*K_const; % Hook's constant.
Mmean = zeros(num_windows, numSweeps_per_split);
Emean = zeros(num_windows, numSweeps_per_split);
X_equilibrium = zeros(1,num_windows);


%% Biased sampling

center_shifts = zeros(num_targets, 1+numUpdates);
estimated_gradient = zeros(1, num_windows);

sim_id = 0;

for target_id = 1:1:num_targets
    sim_id = sim_id + 1;
    centerCurrent = Target_centers(target_id);
    
    X_equilibrium(sim_id) = centerCurrent;
    [Mmean(sim_id,:), Emean(sim_id,:), spin] = MCMC_ising_model(kT, numSweeps_per_split, J, spin, K(target_id), centerCurrent);

    center_shifts(target_id,1) = centerCurrent - nanmean(Mmean(sim_id,:)); % Get the difference between target and center of histogram
    estimated_gradient(sim_id) = 2*K_const*(centerCurrent - nanmean(Mmean(sim_id,:)));
    
    for update_id = 1:(numUpdates)
        sim_id = sim_id + 1;
    
        %%% Update by gradient descent (requires to decide on a learning param)
        [centerUpdated] = update_umbrella_center(Mmean(sim_id-1,:), Target_centers(target_id), centerCurrent, 1); % Used X_eq(i) as Target, and used Mmean(i,:) 
            if centerUpdated < lower_boundary
                centerUpdated = lower_boundary;
            elseif centerUpdated > upper_boundary 
                centerUpdated = upper_boundary;
            end
        %%% Update by analytic solution (via estimated gradient)
        centerUpdated = Target_centers(target_id) + estimated_gradient(sim_id-1)/(2*K_const);
        
        
        centerCurrent = centerUpdated;
        X_equilibrium(sim_id) = centerCurrent;
        [Mmean(sim_id,:), Emean(sim_id,:), spin] = MCMC_ising_model(kT, numSweeps_per_split, J, spin, K(target_id), centerCurrent);
        center_shifts(target_id,1+update_id) = centerCurrent - nanmean(Mmean(sim_id,:)); % Get the difference between target and center of histogram
        estimated_gradient(sim_id) = 2*K_const*(centerCurrent - nanmean(Mmean(sim_id,:)));

        
%         X_updated_record(window_id,1+update_id) = X_equilibrium(window_id);
    end
end

%% Checks

try
    figure(fig_estimatedGradient)
catch
    fig_estimatedGradient = figure();
end
hold on
title("Estimated Gradient from Deviation of histogram centers")
scatter(X_equilibrium,estimated_gradient,'g')
ylabel("Estimated Gradient")
xlabel("RC")

figure; hold on
subplot(2,1,1); hold on
    title("Deviation from the target of sampled histogram centers per update")
    for k = 1:size(center_shifts,2)
        if k == size(center_shifts,2)
            tmp_p = plot(center_shifts(:,k),'go-');
        else
            tmp_p = plot(center_shifts(:,k),'ro-');
            tmp_p.Color(4) = k./size(center_shifts,2);
        end
    end
    xlabel("Umbrella ID"); ylabel("X_{Target} - X_{Histogram}");
subplot(2,1,2); hold on
    for k = 1:size(center_shifts,2)
        if k == size(center_shifts,2)
            tmp_p = plot(abs(center_shifts(:,k)),'go-');
        else
            tmp_p = plot(abs(center_shifts(:,k)),'ro-');
            tmp_p.Color(4) = k./size(center_shifts,2);
        end
    end
    xlabel("Umbrella ID"); ylabel("|X_{Target} - X_{Histogram}|");
    
%%
WHAM_Ising_model()

figure;
scatter(1:length(X_equilibrium),X_equilibrium)