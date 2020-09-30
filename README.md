# Newton-Rhapson-Power-Flow
Implementation of Newton Rhapson Power Flow (NRPF) algorithm on IEEE-14 bus system with transformer taps, Q-limits, and Fast Decoupled power flow approach for better computation. The program is developed without using **any** in-built functions of ```MATLAB```.  

*Note: In order to make this work, you have to assume that Bus 1 is the slack bus.*

![IEEE-14 bus system](./IEEE14bus_data/IEEE14bus.PNG)

# Development Stages
1. Reading bus and branch data in common data format 
1. Y-bus formulation
2. Calculating Jacobian Matrix
3. Inversion using Crout's LU factorization 
4. Solving NRPF unless *power mismatch* < 0.001
5. Implementing Q-lim controls
6. Fast Decoupled Power Flow

# Dependencies
1. ```MATLAB```
Although any version would work, I wrote it in 2020a. Please let me know if the script does not run on your system. I can provide a version compatible script on request.
