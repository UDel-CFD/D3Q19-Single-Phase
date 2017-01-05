# UDel-CFD LBM D3Q19 Repository
This git repository holds Lattice Boltzmann Method D3Q19 Simulation Code. Simulations are designed and tested on NCAR's Yellowstone. Code is compiled with Intel's fortran compiler ifort and Intel's MPI library.

## Current Simulations
* **Channel Flow** - Simulates both turbulent channel and particle-laden channel with MRT collision, swap fused collision and propagation method, and interpolation bounce-back for fluid particle collision.
* **Pipe Flow** - Simulates laminar and turbulent pipe flow with MRT collision, two-array fused collision and propagation method, and interpolation bounce-back for pipe boundary. Branch name: *Pipe-Flow*
* **Urban Canopy** - Simulates turbulence over a urban canopy with BGK collision, two-array fused collision and propagation, and mid-link bounce back on the buildings. Branch name: *Urban-Canopy*

___
## Additional Resources
* [Git Cheat Sheet] [GitCheatSheet]
* [NCAR Yellowstone Website] [NCAR]

   [NCAR]: <https://www2.cisl.ucar.edu/resources/yellowstone>