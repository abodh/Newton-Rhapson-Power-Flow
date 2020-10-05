%{
Newton Rhapson Power Flow
Author: Abodh Poudyal
Last updated: September 30, 2020
%}

clear;
clc;
% format short % to display less significant digits in the result 

%% 1. Reading bus and branch data in common data format
% external function to extract the data from IEEE common data format

bus_path = 'IEEE14bus_data/bus_data.txt';
branch_path = 'IEEE14bus_data/branch_data.txt';
[bus_imp, branch_imp, bus_data, branch_data] = data_extraction(bus_path, branch_path);

%{ 
to reduce the computational complexity, we will only compute 
for existing brances
%}

% define some important variables

% from which bus
from = branch_data(:,1);
% to which bus
to = branch_data(:,2);

% extract voltage data
V_flat = bus_data.data(:,11);  

% flat start means |V| = 1.0 pu and delta = 0
% exact values for PV and slack bus whereas flat start for the rest
V_flat(find(V_flat == 0)) = 1;
delta_flat = zeros(length(V_flat),1);

% V = V_flat;
% delta = delta_flat;

% number of buses in the entire system
n_bus = length(bus_data.data(:,3));

% number of branches
n_branch = length(branch_imp);

% iterate unless power mismatch < 0.01 (tolerance)
tolerance = 0.01;

% base power
base_MW = 100;

%% 2. Calculating the Y-bus matrix

Y_bus = Ybus(n_bus, n_branch, branch_imp, bus_imp, from, to);
G = real(Y_bus); % conductance (G) <- real part of admittance 
B = imag(Y_bus); % susceptance (B) <- the imaginary part of admittance 

%% 3,4,5. Calculating Jacobian Matrix, LU factorization, and NR power flow

Q_lim_status = 1;
while (Q_lim_status)
    % scheduled power
    Ps = (bus_data.data(:,8) - bus_data.data(:,6))/base_MW;
    Qs = (bus_data.data(:,9) - bus_data.data(:,7))/base_MW;

    % number of pq buses
    n_pq = length(find(bus_data.data(:,3) == 0));

    % number of PV buses 
    n_pv = length(find(bus_data.data(:,3) == 2));

    % stores an array of PQ bus IDs
    pq_bus_id = find(bus_data.data(:,3) == 0);

    % stores an array of PV bus IDs 
    pv_bus_id = find(bus_data.data(:,3) == 2);

    % Newton Rhapson Power Flow 
    [Volt, Angle, error_avg] = ...
        NewtonRhapson(tolerance, n_bus, n_pv, n_pq, pq_bus_id, V_flat, ...
        delta_flat, G, B, Ps, Qs);
    
    % Q-limit check
    [Q_lim_status, bus_data] = Qlim(Volt, Angle, bus_data, G, B,...
        base_MW, pv_bus_id, n_bus);
end
    
%% 6. Fast decoupled power flow
[Volt_FD, Angle_FD, error_avg_FD] = ...
    FastDecoupledPF(tolerance, from, to, n_branch, n_bus, n_pv, n_pq, ...
    pq_bus_id, V_flat, delta_flat, G, B, Ps, Qs, branch_imp, bus_imp); 

%% plots of the result

% NRLF Voltage
figure('color', [1,1,1])
subplot(3,2,1)
str = "bus 1";
for i = 1: length(Volt)
    plot(Volt(i,:), 'Linewidth', 1.5)
    hold on
    if i > 1
        str = [str , strcat('bus ' , num2str(i))];
    end
end

ylabel('Voltage (pu)')
xlabel('Number of iteration')
grid on
% set(gca,'XTick',(1:1:4))
set(gca,'gridlinestyle','--','fontname','Times New Roman','fontsize',12);
hold off

% NRLF Angle
subplot(3,2,2)
for i = 1: length(Angle)
    plot(Angle(i,:), 'Linewidth', 1.5)
    hold on
end
ylabel('Angle (rad)')
xlabel('Number of iteration')
grid on
set(gca,'gridlinestyle','--','fontname','Times New Roman','fontsize',12);
hold off

% FDLF Voltage
subplot(3,2,3)
for i = 1: length(Volt_FD)
    plot(Volt_FD(i,:), 'Linewidth', 1.5)
    hold on
end
ylabel('Voltage (pu)')
xlabel('Number of iteration')
grid on
set(gca,'gridlinestyle','--','fontname','Times New Roman','fontsize',12);
hold off

% FDLF Angle
subplot(3,2,4)
for i = 1: length(Angle_FD)
    plot(Angle_FD(i,:), 'Linewidth', 1.5)
    hold on
%     ylim([0:1])
end

ylabel('Angle (rad)')
xlabel('Number of iteration')
grid on
set(gca,'gridlinestyle','--','fontname','Times New Roman','fontsize',12);
hold off

lgd = legend (str, 'NumColumns', 4);
% set(lgd,'position',poshL);      % Adjusting legend's position
% axis(hL,'off');                 % Turning its axis off
