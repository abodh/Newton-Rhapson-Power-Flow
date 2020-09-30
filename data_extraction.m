clear;
clc;
format short % to display less significant digits in the result 

% for convenience, the bus and branch data have been separated to two files
% both of the files are inside IEEE14bus_data folder

%import the data for each buses
bus_data = importdata('IEEE14bus_data/bus_data.txt');

% bus_imp = [conductance susceptance] for each of the bus
bus_imp = bus_data.data(:,14:15); 

%import the data for each branches
branch_data = importdata('IEEE14bus_data/branch_data.txt');

% branch_imp = [R X B tap_setting] for each of the branches
branch_imp = [branch_data(:,7:9), branch_data(:,15)];

%{ 
to reduce the computational complexity, we will only compute 
for existing brances
%}

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

V = V_flat;
delta = delta_flat;

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

% number of branches
n_branch = length(branch_imp);