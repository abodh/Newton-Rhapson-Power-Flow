mismatch_FD = ones(2 * n_bus - n_pv - 2, 1);
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

while any(abs(mismatch_FD(:)) > 0.01)
    iter = iter + 1;
    Volt_FD(:,iter) = V;
    Angle_FD(:,iter) = delta;   
    mismatch_FD = power_mismatch(Ps, Qs, G, B, V, delta, n_bus, pq_bus_id);
    error = croutLU(Jacob_FD, mismatch_FD);
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