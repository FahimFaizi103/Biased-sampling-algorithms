%% Histogram of spin
figure
hold on
for i = 1:1:numSim
histogram(Mmean(i,equilib_time:end))
end
title(['kT = ' num2str(kT)])
xlabel('\langle M \rangle')
%% Probability distribution 
figure
hold on
 for i = 1:1:numSim
 [values,edges]=histcounts(Mmean(i,equilib_time:end),50);
 values = values/sum(values);
 plot(edges(2:length(edges)),values)
 end
title('Probability density')
xlabel('Magnetisation per spin')
ylabel('probability')
title(['kT = ' num2str(kT)])
%% Free energy plot
figure
hold on
Num_Bins = 35;
 for i = 1:1:numSim
     
 [values,edges]=histcounts(Mmean(i,equilib_time:end),Num_Bins);
 values = values/sum(values); 
 plot(edges(2:length(edges)),-log(values)-min(-log(values)),'-','Linewidth',2)
 
 end
 
 title('Free energy profile')
 ylabel('Free energy')
 xlabel('Magnetisation per spin')
 xlim([-1 1])