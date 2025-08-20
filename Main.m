clear; clc;

k = 1;                     % exponential load coefficient
S = [1 10];                % slenderness values
Nx = 500; Ny = 500;
modeList = [0 1 2 50];

nP = numel(S);
for i = 1:nP
    plates(i) = Plate(k, S(i), Nx, Ny);
end

plotter = Plotter(plates, modeList);
plotter.plot_displacement();
plotter.plot_stresses();
