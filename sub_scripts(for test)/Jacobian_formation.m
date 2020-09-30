%{
Author: Abodh Poudyal (11744564)
Washington State University

Last Updated: September 14, 2020 
%}
%%
clear;
clc;
format short % to display less significant digits in the result 

%% data extraction and variable definition from the 14 bus text files
% for convenience, the bus and branch data have been separated to two files

%import the data for each buses
bus_data = importdata('bus_data.txt');

% bus_imp = [conductance susceptance] for each of the bus
bus_imp = bus_data.data(:,14:15); 

%import the data for each branches
branch_data = importdata('branch_data.txt');

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

% number of buses in the entire system
n_bus = length(bus_data.data(:,3));

% number of PV buses 
n_pv = length(find(bus_data.data(:,3) == 2));

% number of pq buses
n_pq = length(find(bus_data.data(:,3) == 0));

% stores an array of PQ bus IDs
pq_bus_id = find(bus_data.data(:,3) == 0);

%% Y-bus calculation

%initialize the Y-bus according to the size of nodes
Y_bus = zeros(length(n_bus));

for i = 1:length(branch_imp)
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

for i  = 1:length(n_bus)
    % the individual buses will have their own shunt device 
    % It should also be included in the Y_bus; Yii = Yii + (Gi + j.Bi)
    Y_bus(i,i) = Y_bus(i,i) + bus_imp(i,1) + 1i*bus_imp(i,2);
end

% end of Y-bus calculation

G = real(Y_bus); % conductance is the real part of admittance
B = imag(Y_bus); % susceptance is the imaginary part of admittance

%% Jacobian 

P = zeros(n_bus,1);
Q = zeros(n_bus,1);

%calculating the active and reactive power at each bus
%{
P(i) = sum(j=1->n) |Vi||Vj|(Gij * cos(delta_i - delta_j) + Bij *
                            sin(delta_i - delta_j)
Q(i) = sum(j=1->n) |Vi||Vj|(Gij * sin(delta_i - delta_j) - Bij *
                            cos(delta_i - delta_j)
%}
for i = 1 : n_bus
    for j = 1 : n_bus
        P(i) = P(i) + V(i)*V(j)*(G(i,j)*cos(delta(i)-delta(j)) + B(i,j)*sin(delta(i)-delta(j)));
        Q(i) = Q(i) + V(i)*V(j)*(G(i,j)*sin(delta(i)-delta(j)) - B(i,j)*cos(delta(i)-delta(j)));
    end
end

%{
slack bus : 1
PV buses = 4
PQ bus = 9
size of Jacobian matrix = 2*n_bus - n_pv - 2 = 22
%}

% size of each of the Jacobian sub matrices
J11 = zeros(n_bus-1, n_bus-1);
J12 = zeros(n_bus-1, n_pq);
J21 = zeros(n_pq, n_bus-1);
J22 = zeros(n_pq, n_pq);

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

Jacob = [J11 J12; J21 J22];

%% calculation of mismatches

Ps = bus_data.data(:,8) - bus_data.data(:,6);
Qs_pseudo = bus_data.data(:,9) - bus_data.data(:,7);
Qs = Qs_pseudo(find(bus_data.data(:,3) == 0));

mismatch = [Ps(2:end);Qs];
V_del = inv(Jacob) * mismatch;


