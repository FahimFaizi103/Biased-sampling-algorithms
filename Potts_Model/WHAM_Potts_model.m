% use this code to unbias the data from biased sampling simulations. 
% This is Weighted Histogram analysis method.
%% Load the data

% load('Potts_12by12.mat')

%% Compute Biased Probability density

Num_Bins = 20;
window = linspace(1,num_windows,num_windows);

figure
hold on
for i = 1:1:num_windows
    
    [values,edges]=histcounts(Mmean(window(i),:),Num_Bins);
    values = values/sum(values);
    
    plot(edges(2:end), values, 'Linewidth', 2) % probability distribution
%     plot(edges(2:end), -log(values)-min(-log(values)), 'Linewidth',2) % free energy
    
end

box on
xlabel('$m_0$','Interpreter','latex')
ylabel('probability')
xlim([min_q, max_q])
%% Computing F_i values

beta = 1./kT;
diff_threshold = 10^(-3);
diff = 10*diff_threshold;

F_record = [];
F = ones(1,num_windows);
F_record(1,:) = F;
%%
tt = 1;
while diff > diff_threshold
    
    tt = tt+1;
    
    bias = zeros(num_windows, numSweeps);
    numerator = zeros(num_windows, numSweeps, num_windows);
    denominator = zeros(num_windows, numSweeps);
    
    for kk = 1:1:num_windows
        
        bias = compute_bias(K(kk), Mmean, X_equilibrium(kk));
        numerator(:,:,kk) = exp(-beta.*bias);
        
        denominator = denominator + numSweeps.*exp(-beta.*bias + beta.*F(kk));
        
    end
    
    arg = numerator./denominator;
    arg = sum(sum(arg,2));
    
    F = -kT.*log(arg(:));
    F_record = [F_record; F'];
    
    diff = abs(max(F_record(tt,:) - F_record(tt-1,:)));
    
    
end

%% plot the F_i values against num_iterations to visualize convergence

 figure
 plot(F_record)
 ylabel('$F_i$','interpreter','latex')
 xlabel('iteration')
 %% Computing the unbiased probability distribution.
 
% load('umbrella_sampling_1.mat')
% load('data','F')

Num_Bins = 25;
EDGES = linspace(lower_boundary,upper_boundary,Num_Bins+1);
beta = 1./kT;
prob = zeros();
for i = 2:1:length(EDGES)

    if i == length(EDGES)
        
    candidates = find(Mmean >= EDGES(i-1) & Mmean <= EDGES(i) );
    
    else
        
    candidates = find(Mmean >= EDGES(i-1) & Mmean < EDGES(i) );
        
    end
    
    Sum = 0;
    for kk = 1:1:length(candidates)
        
        denom = 0;
        for j = 1:1:num_windows
            
            bias = compute_bias(K(j), EDGES(i), X_equilibrium(j));
            denom = denom + numSweeps.*exp(-beta.*bias + beta.*F(j));
                      
        end
        
        Sum = Sum + 1./denom;
              
    end
    
    prob(i-1) = Sum;
    
end

prob = prob./sum(prob);

%% Plot the probability distribution.

figure
hold on
plot(EDGES(2:end), prob, 'Linewidth', 2) % probability distribution
xlabel('$m$','interpreter','latex')
ylabel('probability')
%% Plot the free energy profile

figure
hold on
plot(EDGES(2:end), -log(prob)-min(-log(prob)), 'Linewidth',2)
box on
xlabel('$m$','interpreter','latex')
ylabel('Free energy')
