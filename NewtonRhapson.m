% Power flow iteration
iter = 0;
mismatch_NR = ones(2 * n_bus - n_pv - 2, 1);

while any(abs(mismatch_NR(:)) > 0.01)
    iter = iter + 1;
    Volt(:,iter) = V;
    Angle(:,iter) = delta;
    Jacob = Jacobian(V, delta, n_bus, n_pq, pq_bus_id, G, B);
    mismatch_NR = power_mismatch(Ps, Qs, G, B, V, delta, n_bus, pq_bus_id);
%   error = inv(Jacob) * mismatch;
    error = croutLU(Jacob, mismatch_NR);
    delta(2:end) = delta(2:end) + error(1 : n_bus-1);
    V(pq_bus_id) = V(pq_bus_id) + error(n_bus : end);
end

function mismatch_NR = power_mismatch(Ps, Qs, G, B, V, delta, n_bus, pq_bus_id)
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
    mismatch_NR = [delta_P(2:end);delta_Q];    
end