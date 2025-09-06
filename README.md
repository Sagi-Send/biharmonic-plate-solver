# Biharmonic Plate Solution

This repository implements a numerical solver for the in-plane displacements and stresses of a propped-cantilever rectangular plate subjected to a transverse continuous load. The formulation integrates the plane-stress constitutive laws, enforces shear compatibility, and applies weak clamping.

## Theory

For a linearly elastic plate, the Airy stress function $\Phi(x,y)$ satisfies the biharmonic equation
$$
\nabla^{4}\,\Phi(x,y) \;=\; q(x,y).
$$
Here the transverse load is taken as an exponential distribution along the span,
$$
w(x) \;=\; w_0\,e^{k x/\ell}.
$$
Stresses are recovered from $\Phi$, and in-plane displacements $u(x,y)$ and $v(x,y)$ follow from integrating the plane-stress strains. The series solution is truncated to a specified number of Fourier modes (see [1,2]).

## Running the code

The scripts are written for MATLAB. To reproduce figures:

```matlab
% at the MATLAB prompt
>> Main
