# Biharmonic Plate Solution

This repository implements a numerical solver for the in-plane displacements and stresses of a propped-cantilever rectangular plate subjected to a transverse continuous load. The formulation integrates the plane-stress constitutive laws, enforces shear compatibility, and applies weak clamping.

## Theory

For a linearly elastic plate, the Airy stress function $\Phi(x,y)$ satisfies the biharmonic equation

$$\nabla^{4}\Phi(x,y)=q(x,y).$$

Here the transverse load is taken as an exponential distribution along the span,

$$w(x) = w_0e^{k x/\ell}.$$

Stresses are recovered from $\Phi$, and in-plane displacements $u(x,y)$ and $v(x,y)$ follow from integrating the plane-stress strains.


## Running the code

The scripts are written for MATLAB.  To reproduce figures:

```matlab
% at the MATLAB prompt
>> Main
```

`Main.m` constructs `Plate` objects for several slenderness values and calls the plotting routines.  Additional plots such as `plot_displacement_exp` or `plot_stresses_exp` can be enabled in `Main.m`.

## Example output

![Deformation field](Figures/cont_deformation.png)

![Superposed stresses](Figures/stress-superpose.png)

## References

1. S. Timoshenko and S. Woinowsky-Krieger, *Theory of Plates and Shells*, McGraw–Hill, 1959.
2. J. N. Reddy, *Theory and Analysis of Elastic Plates and Shells*, CRC Press, 2006.
