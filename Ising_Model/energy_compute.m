function[Emean] = energy_compute(spin,J)


     sumOfNeighbors = ...
      circshift(spin, [ 0  1]) ...
    + circshift(spin, [ 0 -1]) ...
    + circshift(spin, [ 1  0]) ...
    + circshift(spin, [-1  0]);

% each spin sigma_m contributes an energy Em to the total energy.
Em = -J*spin.*sumOfNeighbors; 

% The total enery is of course the sum over each dimension of Em

E = 0.5*sum(sum(Em));
Emean = E./numel(spin); % energy per spin

end