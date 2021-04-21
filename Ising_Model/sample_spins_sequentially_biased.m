function[dE, row, col] = sample_spins_sequentially_biased(spin,J,i,K, X_equilibrium)

 % Picking a sample where sample = ith linear index of matrix: spin.
sample = i;

% LinearIndex picks a spin by acing as an index into the Matrix: spin. We
% can convert it into [row, col] of our matrix spin by the following:

[row,col] = ind2sub(size(spin),sample);

% The subscripts of the neighbours of spin(row,col) are the following.
% The "above" and "below" identify the rows that the nearest neighbours of
%spin(row,col) lies, likewise the "left" and "right" identify the coloums in
%which the nearest neighbours lie.

above = mod(row - 1 - 1, size(spin, 1)) + 1;
below = mod(row + 1 - 1, size(spin, 1)) + 1;
left  = mod(col - 1 - 1, size(spin, 2)) + 1;
right = mod(col + 1 - 1, size(spin, 2)) + 1;

% We now specify the coordinates of the nearest neighbours. 
neighbors = [       spin(above, col);
    spin(row, left);                 spin(row, right);
                    spin(below, col)];

U_bf = K.*(mean(spin(:)) - (2.*spin(row,col))./numel(spin) - X_equilibrium).^2;
U_bi = K.*(mean(spin(:)) - X_equilibrium).^2;

% Change in energy if this spin flips
dE = 2*J.*spin(row,col).*sum(neighbors) + (U_bf - U_bi);

end