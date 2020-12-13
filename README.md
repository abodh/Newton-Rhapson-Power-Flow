# Newton-Rhapson-Power-Flow
Implementation of Newton Rhapson Power Flow (NRPF) algorithm on IEEE-14 bus system with transformer taps, Q-limits, and Fast Decoupled power flow approach for better computation. The program is developed without using **any** in-built functions of ```MATLAB```. Although presented for the IEEE-14 bus system, the code is generalized to work for different system provided that the data is presented in IEEE common data format. Before running it on a test case, please separate bus and branch data in a separate text file. The program evaluates the Q-limit for Newton-Rhapson but not for fast decoupled power-flow. The Q-limit function can easily be integrated with fast decoupled power-flow algorithm. You just need to call the function after being done with the power flow.  

![IEEE-14 bus system](./IEEE14bus_data/IEEE14bus.PNG)

# Assumptions
1. In order to make this work, we have to assume that Bus 1 is the slack bus. Hence, bus 1 is always considered as a slack bus.
2. The Crout's LU factorization algorithm does not work when the first diagonal element of L is zero.

# Development Stages
1. Reading bus and branch data in common data format 
2. Y-bus formulation
3. Calculating Jacobian Matrix
4. Crout's LU factorization to solve unknown values 
5. Solving NRPF unless *power mismatch* < 0.01
6. Fast Decoupled Power Flow
  * when the B matrices are same for both the diagonal sub-Jacobian matrices
  * when the B matrices are different for the diagonal sub-Jacobian matrices i.e. when lossless system is assumed (i.e. B = 1/X).
7. Implementing Q-lim controls
  * evaluated once the power-flow is completed
  * does not evaluate for each iteration

# Dependencies
1. ```MATLAB```
Although any version would work, I wrote it in 2020a. Please let me know if the script does not run on your system. I can provide a version compatible script on request.

# Update
This tool now also solves an IEEE-30 bus system. For higher number of buses there are more than 1 slack buses. Hence, with more than one slack buses, this tool does not provide any power flow solution.
