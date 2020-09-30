%{
Newton Rhapson Power Flow
Author: Abodh Poudyal
Last updated: September 30, 2020
%}

clear;
clc;
format short % to display less significant digits in the result 

% for convenience, the bus and branch data have been separated to two files

%import the data for each buses
bus_data = importdata('IEEE14bus_data/bus_data.txt');

% bus_imp = [conductance susceptance] for each of the bus
bus_imp = bus_data.data(:,14:15); 

%import the data for each branches
branch_data = importdata('IEEE14bus_data/branch_data.txt');

% bus_imp = [R X B tap_setting] for each of the branches
branch_imp = [branch_data(:,7:9), branch_data(:,15)];

% to reduce the computational complexity, we will only consider existing branches
from = branch_data(:,1);
to = branch_data(:,2);

% extracts the voltage data
% exact for PV and slack and flat start for the rest
V = bus_data.data(:,11);  

% flat start means |V| = 1.0 pu and delta = 0
V(find(V == 0)) = 1;
delta = zeros(length(V),1);

V_flat = V;
delta_flat = delta;

% scheduled power
Ps = (bus_data.data(:,8) - bus_data.data(:,6))*0.01;
Qs = bus_data.data(:,9) - bus_data.data(:,7)*0.01;

% number of buses in the entire system
n_bus = length(bus_data.data(:,3));

% number of PV buses 
n_pv = length(find(bus_data.data(:,3) == 2));

% number of pq buses
n_pq = length(find(bus_data.data(:,3) == 0));

% stores an array of PQ bus IDs
pq_bus_id = find(bus_data.data(:,3) == 0);

n_branch = length(branch_imp);

Y_bus = Ybus(n_bus, n_branch, branch_imp, bus_imp, from, to);
G = real(Y_bus); % conductance is the real part of admittance
B = imag(Y_bus); % susceptance is the imaginary part of admittance
iter = 0;

% Power flow iteration
mismatch = ones(2 * n_bus - n_pv - 2, 1);

%% Fast decoupled power flow
mismatch = ones(2 * n_bus - n_pv - 2, 1);
iter = 0;
V = V_flat;
delta = delta_flat;
% B11 = B(2:end,2:end);

B_prime = zeros(n_bus,n_bus);
for i = 1 : n_branch
    % Yij = -((1/(rij + j.xij)) + j.Bij/2)/tap-setting 
    B_prime(from(i),to(i)) = - (1/(0 + 1i*(branch_imp(i,2))))/branch_imp(i,4);

    %Yij = Yji
    B_prime(to(i),from(i)) = B_prime(from(i),to(i));

    % considering tap-setting the admittance matrix looks like
    % Y = [Y/t^2  -Y/t
    %       -Y/t     Y]
    B_prime(from(i),from(i)) = B_prime(from(i),from(i)) + ((1/(0 + 1i*(branch_imp(i,2))))/(branch_imp(i,4))^2) + 1i*0.5*branch_imp(i,3);

    B_prime(to(i),to(i)) = B_prime(to(i),to(i)) + (1/(0 + 1i*(branch_imp(i,2)))) + 1i*0.5*branch_imp(i,3);
end

for i  = 1 : n_bus
    % the individual buses will have their own shunt device 
    % It should also be included in the Y_bus; Yii = Yii + (Gi + j.Bi)
    B_prime(i,i) = B_prime(i,i) + 0 + 1i*bus_imp(i,2);
end
B_prime_ = imag(B_prime);
B11 = B_prime_(2:end,2:end);

for i = 1 : length(pq_bus_id)
    for j = 1 : length(pq_bus_id)
        B22(i,j) = B(pq_bus_id(i), pq_bus_id(j));
    end
end

Jacob_FD = -[B11 zeros(n_bus-1, length(pq_bus_id)); zeros(length(pq_bus_id),n_bus-1) B22];

while any(abs(mismatch(:)) > 0.01)
    iter = iter + 1;
    Volt_FD(:,iter) = V;
    Angle_FD(:,iter) = delta;   
    mismatch = power_mismatch(Ps, Qs, G, B, V, delta, n_bus, pq_bus_id);
    error = croutLU(Jacob_FD, mismatch);
    delta(2:end) = delta(2:end) + error(1 : n_bus-1);
    V(pq_bus_id) = V(pq_bus_id) + error(n_bus : end);
end
for i = 1: length(Volt_FD)
    plot(Volt_FD(i,:))
    hold on
%     ylim([0:1])
end
hold off

figure
for i = 1: length(Angle_FD)
    plot(Angle_FD(i,:))
    hold on
%     ylim([0:1])
end
%%

while any(abs(mismatch(:)) > 0.01)
    iter = iter + 1;
    Volt(:,iter) = V;
    Angle(:,iter) = delta;
    Jacob = Jacobian(V, delta, n_bus, n_pq, pq_bus_id, G, B);
    mismatch = power_mismatch(Ps, Qs, G, B, V, delta, n_bus, pq_bus_id);
%   error = inv(Jacob) * mismatch;
    error = croutLU(Jacob, mismatch);
    delta(2:end) = delta(2:end) + error(1 : n_bus-1);
    V(pq_bus_id) = V(pq_bus_id) + error(n_bus : end);
end

% plot (avg_error,'r', 'Linewidth',2)
% ylabel('Averaged absolute error')
% xlabel('number of iteration')
% grid on
% set(gca,'XTick',(1:1:4))
% set(gca,'gridlinestyle','--','fontname','Arial','fontsize',12);

for i = 1: length(Volt)
    plot(Volt(i,:))
    hold on
%     ylim([0:1])
end
hold off

figure
for i = 1: length(Volt)
    plot(Angle(i,:))
    hold on
%     ylim([0:1])
end

function error_VD = croutLU(A, b)
L = zeros(length(A));
U = zeros(length(A));

L(1:length(A),1) = A(1:length(A),1);
for i = 1:length(A)
    U(i,i) = 1;
end
U(1,2:length(A)) = A(1,2:length(A))/L(1,1);

for j = 2 : length(A)
    for k = j : length(A)
        L(k,j) = A(k,j) - L(k,1:j-1)*U(1:j-1,j);
    end
    for k = j+1 : length(A)
        U(j,k) = (A(j,k) - (L(j,1:j-1)*U(1:j-1,k)))/L(j,j);
    end
end

% forward substitution
y = zeros(length(A),1);
for i = 1:length(A)
        if i == 1
            y(1) = b(1)/L(1,1);
        else
            y(i) = (b(i)-L(i,1:i-1)*y(1:i-1))/L(i,i);
        end
end

% backward substitution
x = zeros(length(A),1);
for i = length(A):-1:1   
    if i == length(A)
        x(length(A)) = y(length(A))/U(length(A),length(A));
    else
        x(i) = (y(i)-U(i,i+1:length(A))*x(i+1:length(A)))/U(i,i);
    end
end
error_VD = x;
end


function mismatch = power_mismatch(Ps, Qs, G, B, V, delta, n_bus, pq_bus_id)
    P = zeros(n_bus,1);
    Q = zeros(n_bus,1);
    
    % calculating the active and reactive power at each bus
    %{
        P(i) = sum(j=1->n) |Vi||Vj|(Gij * cos(delta_i - delta_j) + 
                                    Bij * sin(delta_i - delta_j)
        Q(i) = sum(j=1->n) |Vi||Vj|(Gij * sin(delta_i - delta_j) - 
                                    Bij * cos(delta_i - delta_j)
    %}
    for i = 1 : n_bus
        for j = 1 : n_bus
            P(i) = P(i) + V(i)*V(j)*(G(i,j)*cos(delta(i)-delta(j)) + B(i,j)*sin(delta(i)-delta(j)));
            Q(i) = Q(i) + V(i)*V(j)*(G(i,j)*sin(delta(i)-delta(j)) - B(i,j)*cos(delta(i)-delta(j)));
        end
    end
    
    delta_P = Ps - P;
    delta_Q = Qs - Q;
    delta_Q = delta_Q(pq_bus_id);
    mismatch = [delta_P(2:end);delta_Q];    
end


function Y_bus = Ybus(n_bus, n_branch, branch_imp, bus_imp, from, to)
    % Y-bus calculation
    % initialize the Y-bus according to the size of nodes
    Y_bus = zeros(n_bus);

    for i = 1 : n_branch
        % Yij = -((1/(rij + j.xij)) + j.Bij/2)/tap-setting 
        Y_bus(from(i),to(i)) = - (1/(branch_imp(i,1) + 1i*(branch_imp(i,2))))/branch_imp(i,4);

        %Yij = Yji
        Y_bus(to(i),from(i)) = Y_bus(from(i),to(i));

        % considering tap-setting the admittance matrix looks like
        % Y = [Y/t^2  -Y/t
        %       -Y/t     Y]
        Y_bus(from(i),from(i)) = Y_bus(from(i),from(i)) + ((1/(branch_imp(i,1) + 1i*(branch_imp(i,2))))/(branch_imp(i,4))^2) + 1i*0.5*branch_imp(i,3);
                                 
        Y_bus(to(i),to(i)) = Y_bus(to(i),to(i)) + (1/(branch_imp(i,1) + 1i*(branch_imp(i,2)))) + 1i*0.5*branch_imp(i,3);
    end

    for i  = 1 : n_bus
        % the individual buses will have their own shunt device 
        % It should also be included in the Y_bus; Yii = Yii + (Gi + j.Bi)
        Y_bus(i,i) = Y_bus(i,i) + bus_imp(i,1) + 1i*bus_imp(i,2);
    end
end

function Jacob_matrix = Jacobian(V, delta, n_bus, n_pq, pq_bus_id, G, B)
    P = zeros(n_bus,1);
    Q = zeros(n_bus,1);

    % calculating the active and reactive power at each bus
    %{
        P(i) = sum(j=1->n) |Vi||Vj|(Gij * cos(delta_i - delta_j) + 
                                    Bij * sin(delta_i - delta_j)
        Q(i) = sum(j=1->n) |Vi||Vj|(Gij * sin(delta_i - delta_j) - 
                                    Bij * cos(delta_i - delta_j)
    %}
    for i = 1 : n_bus
        for j = 1 : n_bus
            P(i) = P(i) + V(i)*V(j)*(G(i,j)*cos(delta(i)-delta(j)) + B(i,j)*sin(delta(i)-delta(j)));
            Q(i) = Q(i) + V(i)*V(j)*(G(i,j)*sin(delta(i)-delta(j)) - B(i,j)*cos(delta(i)-delta(j)));
        end
    end

    %{
    For IEEE 14 bus system
    slack bus : 1
    PV buses = 4
    PQ bus = 9
    size of Jacobian matrix = 2*n_bus - n_pv (4) - 2 = 22
    %}

    % size of each of the Jacobian sub matrices
    J11 = zeros(n_bus-1, n_bus-1);
    J12 = zeros(n_bus-1, n_pq); % remove the bus (cols) where V is given
    J21 = zeros(n_pq, n_bus-1); % remove the bus (rows) where Q is unknown
    J22 = zeros(n_pq, n_pq);
    
    %{
    refer to the link below for sub-matrices equations:
    https://bit.ly/2GU1Xx0
    
    B. Sereeter, C. Vuik, and C. Witteveen
    REPORT 17-07 
    On a comparison of Newton-Raphson solvers for power flow problems
    %}
    
    % Calculating J11
    for i = 2 : n_bus
        for k = 2 : n_bus
            if (i == k)
                J11(i-1,k-1) = - Q(i) - (V(i)^2 * B(i,i));
            else
                J11(i-1,k-1) = V(i)*V(k)*(G(i,k)*sin(delta(i)-delta(k)) - B(i,k)*cos(delta(i)-delta(k)));
            end
        end
    end

    % Calculating J21
    for i = 2 : n_pq + 1 
        j = pq_bus_id(i-1);
        for k = 2 : n_bus
            if (j == k)
                J21(i-1,k-1) = P(j) - (V(j)^2 * G(j,j));
            else
                J21(i-1,k-1) = -V(j)*V(k)*(G(j,k)*cos(delta(j)-delta(k)) + B(j,k)*sin(delta(j)-delta(k)));
            end
        end
    end


    % Calculating J12
    for i = 2 : n_bus  
        for k = 2 : n_pq + 1
            j = pq_bus_id(k-1);
            if (i == j)
                J12(i-1,k-1) = P(j) + (V(j)^2 * G(j,j));
            else
                J12(i-1,k-1) = V(i)*(G(i,j)*cos(delta(i)-delta(j)) + B(i,j)*sin(delta(i)-delta(j)));
            end
        end
    end

    % Calculating J22
    for i = 2 : n_pq + 1
        j = pq_bus_id(i-1);
        for k = 2 : n_pq + 1
            l = pq_bus_id(k-1);
            if (j == l)
                J22(i-1,k-1) = Q(j) - (V(j)^2 * B(j,j));
            else
                J22(i-1,k-1) = V(j)*(G(j,l)*sin(delta(j)-delta(l)) - B(j,l)*cos(delta(j)-delta(l)));
            end
        end
    end
    
    % combining the jacobian matrix
    Jacob_matrix = [J11 J12; J21 J22];    
end


