% [X_updated] = update_umbrella_center(Mmean(i,:), X_equilibrium(i)); % Used X_eq(i) as Target, and used Mmean(i,:) 

function [centerUpdated] = update_umbrella_center(positionsSampled, targetCenter, currentCenter, learningParam)
    centerShift = targetCenter - nanmean(positionsSampled); % Get the difference between target and center of histogram
    centerUpdated = currentCenter + learningParam*centerShift;
end