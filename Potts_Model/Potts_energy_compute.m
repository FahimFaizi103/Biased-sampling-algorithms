function[Emean] = Potts_energy_compute(spin,J)

E_m = zeros();


for i = 1:1:numel(spin)
    
    spin_sample = i;
    
    [row,col] = ind2sub(size(spin),spin_sample);
    
    % The subscripts of the neighbours of spin(row,col) are the following.
    % The "above" and "below" identify the rows that the nearest neighbours of
    %spin(row,col) lies, likewise the "left" and "right" identify the coloums in
    %which the nearest neighbours lie.
    
    above = mod(row - 1 - 1, size(spin, 1)) + 1;
    below = mod(row + 1 - 1, size(spin, 1)) + 1;
    left  = mod(col - 1 - 1, size(spin, 2)) + 1;
    right = mod(col + 1 - 1, size(spin, 2)) + 1;

    % We now specify the state of the nearest neighbours. 
     neighbors = [       spin(above, col);
         spin(row, left);                 spin(row, right);
                         spin(below, col)];
    
    neighbors(neighbors ~= spin(row,col)) = 0; % unlike neighbours contribute zeros to Hamiltonian
    neighbors(neighbors == spin(row,col)) = 1;
    
    E_m(i) = -J.*sum(neighbors);
    
    
end

E = 0.5.*sum(E_m); % total energy.
Emean = E./numel(spin);

end
