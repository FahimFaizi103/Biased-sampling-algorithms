
clear all
clc

%% load the data

load('12by12_ising_kT_2_standard_MCMC.mat')

%% Histogram of spin

figure
hold on  

histogram(Mmean(equilib_time:end))

title(['kT = ' num2str(kT)])
xlabel('m')
ylabel('Frequency')
box on
%% Probability distribution

figure
hold on
Num_Bins = 48;

[values,edges]=histcounts(Mmean(equilib_time:end),Num_Bins);
values = values/sum(values);
plot(edges(2:length(edges)),values,'-','Linewidth',2)


title('Probability density')
xlabel('m') % magnetisation density i.e. m = M/N.
ylabel('probability')
title(['kT = ' num2str(kT)])
xlim([-1 1])
box on
%% Free energy plot

figure
hold on
Num_Bins = 48;

[values,edges]=histcounts(Mmean(equilib_time:end),Num_Bins);
values = values/sum(values);
plot(edges(2:length(edges)),-log(values)-min(-log(values)),'-','Linewidth',2)

title('Free energy profile')
ylabel('\beta F = -ln(Z)')
xlabel('m')
xlim([-1 1])
box on