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

%% Y-bus calculation

%initialize the Y-bus according to the size of nodes
Y_bus = zeros(length(bus_imp));

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

for i  = 1:length(bus_imp)
    % the individual buses will have their own shunt device 
    % It should also be included in the Y_bus; Yii = Yii + (Gi + j.Bi)
    Y_bus(i,i) = Y_bus(i,i) + bus_imp(i,1) + 1i*bus_imp(i,2);
end
csvwrite('Y_bus_abodh.csv', Y_bus) 

%% Some additional tests (ignore)
% Y_inter = 0;
% c = 1;
%     while (c < 21 && i==from(c))
%         Y_inter = Y_inter + ((1/(branch_imp(c,1) + 1i*(branch_imp(c,2))) + 1i*0.5*branch_imp(c,3)));
%         c = c + 1;
%     end
%     Y_bus(i,i) = Y_inter;
%     Y_inter = 0;

