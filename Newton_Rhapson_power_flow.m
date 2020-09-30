%{
Newton Rhapson Power Flow
Author: Abodh Poudyal
Last updated: September 30, 2020
%}

% 1. Reading bus and branch data in common data format
% external function to extract the data from IEEE common data format
data_extraction();

% 2. Calculating the Y-bus matrix
Y_bus = Ybus(n_bus, n_branch, branch_imp, bus_imp, from, to);
G = real(Y_bus); % conductance (G) <- real part of admittance 
B = imag(Y_bus); % susceptance (B) <- the imaginary part of admittance 

% 3. Calculating Jacobian Matrix
% 4. Crout's LU factorization to solve unknown values
% 5. Solving NRPF unless power mismatch < 0.001
NewtonRhapson();

% 6. Fast decoupled power flow
FastDecoupledPF();

% 7. Analyzing the Q-limits Q_min < Q < Q_max

% plots of the result
for i = 1: length(Volt)
    plot(Volt(i,:), 'Linewidth', 2)
    hold on
%     ylim([0:1])
end
ylabel('Number of iteration')
xlabel('Voltage (pu)')
grid on
set(gca,'XTick',(1:1:4))
set(gca,'gridlinestyle','--','fontname','Times New Roman','fontsize',12);
hold off

figure
for i = 1: length(Volt)
    plot(Angle(i,:))
    hold on
%     ylim([0:1])
end














