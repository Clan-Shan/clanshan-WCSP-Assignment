function w = msinr(a_theta0, x, sp, a, Thetas)
    [M, K] = size(x);

    s_hat = a_theta0 * sp;  % M x K

    % Covariance matrices
    R_in = (x - s_hat)*(x - s_hat)' / K;   % Rn (interference + noise)
    Rs = s_hat * s_hat' / K;               

    % Regularization (to avoid singular matrix)
    epsilon = 1e-6;
    R_in = R_in + epsilon * eye(M);

    % Solve generalized eigenvalue problem
    [V, D] = eig(Rs, R_in);
    [~, idx] = max(real(diag(D)));
    w_temp = V(:, idx);

    % Phase alignment (so 0° will be peak)
    [~, idx0] = min(abs(Thetas - 0));  % 找最接近 0° 的 index
    ref = a(:, idx0);
    phase_ref = angle(w_temp' * ref);
    w = w_temp * exp(-1j * phase_ref);  % Return aligned beamformer
end

