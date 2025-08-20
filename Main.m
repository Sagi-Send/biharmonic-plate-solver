clear; clc;

k = 1;                 % exponential load exponent
S = 1;                 % slenderness S = l/(2h)  (set to 1 ⇒ l = 2h)
Nx = 500; Ny = 500;

plate    = Plate(k, S, Nx, Ny);
modeList = [0 1 2 50];

plotter = Plotter(plate, modeList);
plotter.plot_displacement();
plotter.plot_stresses();
