function [bus_imp, branch_imp, bus_data, branch_data] = ...
    data_extraction (bus_path, branch_path)

    % the bus and branch data have been separated to two files
    % both of the files are inside IEEE14bus_data folder

    %import the data for each buses
    bus_data = importdata(bus_path);

    % bus_imp = [conductance susceptance] for each of the bus
    bus_imp = bus_data.data(:,14:15); 

    %import the data for each branches
    branch_data = importdata(branch_path);

    % branch_imp = [R X B tap_setting] for each of the branches
    branch_imp = [branch_data(:,7:9), branch_data(:,15)];
end