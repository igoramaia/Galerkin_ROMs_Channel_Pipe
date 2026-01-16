# Galerkin ROMs for turbulent channel and pipe flows

Galerkin_ROMs_Channel_Pipe contains a set of codes to generate, run and explore Galerkin-based reduced-order models (ROMs) for turbulent channel and pipe flows. The models use controllability or balanced modes as basis functions, which are derived from the response of the linearised Navier-Stokes system to white-noise forcing. In our associated repository "Galerkin_ROM_CouetteFlow", ROMs for turbulent Couette flow are derived based on controllability modes, with the controllability Gramian being obtained from a Lyapunov equation, assuming a system in statistical steady state, and the controllability modes computed from its eigendecomposition. The channel and pipe flow ROMs described here adopt a different strategy for basis obtention: controllability and observability Gramians are approximated using resolvent modes and a snapshot-based technique, following the methodology laid out in [1,2]. From that approximation, controllability and balanced modes are computed. Those modes are used to perform Galerkin projections and run a system of ODEs for the basis expansion coefficients, as described in [3], allowing subsequent reconstruction of the velocity fields.


ROM generation, time integration and post-processing are performed in the Matlab environment. The codes use a Fourier-Chebyshev discretisation of the domain. Time integration is carried out either with and a standard 4th/5th Runge-Kutta method or with a variable-order solver. The codes can be downloaded from the terminal by typing:

> `git clone https://github.com/igoramaia/Galerkin_ROMs_Channel_Pipe.git`

The files in the root of each folder are the main scripts to generate the ROMs (pre-computation of the operators) and run the system of ODE's. Auxiliary functions to perform mathematical operations are located in the **aux_functions** folder.

# References <h3>

1. G. Dergham, D. Sipp, and J-Ch. Robinet. "Stochastic dynamics and model reduction of amplifier flows: the backward facing step flow". Journal of Fluid Mechanics, 719:406-430, 2013. https://doi.org/10.1017/jfm.2012.610
2. A. V. G. Cavalieri. "Non-linear Galerkin reduced-order models of a mixing layer". AIAA AVIATION 2023 Forum. https://doi.org/10.2514/6.2023-4483
3. A. V. G. Cavalieri and P. A. S. Nogueira. "Reduced-order Galerkin models of plane Couette flow". Physical Review Fluids, 7, L102601, 2022. https://doi.org/10.1103/PhysRevFluids.7.L102601

