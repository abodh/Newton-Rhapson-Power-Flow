# Newton-Rhapson-Power-Flow
Implementation of Newton Rhapson Power Flow (NRPF) algorithm on IEEE-14 bus system with transformer taps, Q-limits, and fast decoupled approach for better computation. The program is developed without using any in-built functions of MATLAB, including inversion.  

# Development Stages
1. Reading bus and branch data in common data format 
1. Y-bus formulation
2. Calculating Jacobian Matrix
3. Inversion using Crout's LU factorization 
4. Solving NRPF unless error < 0.001
5. Implementing Q-lim controls
6. Fast Decoupled Power Flow
