function[spin] = spin_flip_Met(dE, kT, spin, row, col)


P = exp(-dE./kT); % Metropolis criteria
u = rand(1);
% if dE < 0 accept the spin flip, else accept it with probability P
if u <=  min(1,P)
    spin(row,col) = -spin(row,col);
else
    spin(row,col) = spin(row,col);  
end


end