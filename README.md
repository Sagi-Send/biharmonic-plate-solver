# Biharmonic Plate Solution

This repository implements a numerical solver for the in-plane displacements and stresses of a propped-cantilever rectangular plate subjected to a transverse continuous load.

![Problem Formulation](Figures/Problem Formulation.png)

## Theory

For a linearly elastic plate, the Airy stress function $\Phi(x,y)$ identically satisfies equillibrium and enforces compatibiltiy via the biharmonic equation

$$\nabla^{4}\Phi(x,y)=q(x,y).$$

Stresses are recovered from $\Phi$. In-plane displacements $u(x,y)$ and $v(x,y)$ follow from integrating the plane-stress strains and enforcing shear compatibility.
Clamping is applied weakly via either $\partial v(x,y)/\partial x$ or alternatively $\partial u(x,y)/\partial y$. Thus, clamping is applied either at the horizontal or vertical filament, respectively.

![Clamping](Figures/clamping conditions.PNG)

## Resulting Deformation Fields
\textbf{Cantilever}
![Deformation field](Figures/cont_deformation.png)

\textbf{Propped Cantilever}
![Superposed stresses](Figures/stress-superpose.png)

## Running the code

The scripts are written for MATLAB.  To reproduce figures:

```matlab
% at the MATLAB prompt
>> Main
```

`Main.m` constructs `Plate` objects for several slenderness values and calls the plotting routines.  Additional plots such as `plot_displacement_exp` or `plot_stresses_exp` can be enabled in `Main.m`.

## References

1. S. Timoshenko and S. Woinowsky-Krieger, *Theory of Plates and Shells*, McGraw–Hill, 1959.
2. J. N. Reddy, *Theory and Analysis of Elastic Plates and Shells*, CRC Press, 2006.
