# Adaptive-Mesh-Refinement

This directory contains all the code required to perform adaptive mesh refinement using the AMReX software on the simulations of Euler equations. Final report `AMREX_report.pdf` submitted has also been added.

AMReX - https://amrex-codes.github.io/amrex/

* **convergence_amr.cpp** - This file contains the code for performing the convergence analysis on the simulation results. 
* **RiemannExactSol.cpp** - This file contains the solver for computing the exact results of 1-D Euler.
* **eos_v2.cpp** - This file contains the ideal gas equation of state functions. 
* **eos_v2.H** -  This file contains the function headers for `eos_Euler.cpp`.
* **AmrLevelAdv.cpp** - This originated 1-D advection code from Dr. Stephen Millmore has been modified to work on computing 2-D Euler equations. Changes can be tracked by following code labelled with '2020H'.
* **AmrLevelAdv.H** - This file contains the function headers for `AmrLevelAdv.cpp`.
* **hllc.cpp** - This file contains the functions for SLIC numerical solver.  
* **hllc.H** - This file contains the function headers for `hllc.cpp`.
* **Tagging_nd.f90** - This file contains the FORTRAN code for tagging cells for refinement based on gradient calculation.
~`
