function rhs = Kuramoto_RHS(psi,omega,K)
    N = length(psi);
    rhs = zeros(size(psi));
    for j = 1:N
        rhs(j) = omega(j);
        for k = 1:N
            rhs(j) = rhs(j) + (K/N)*sin(psi(k)-psi(j));
        end
    end
end