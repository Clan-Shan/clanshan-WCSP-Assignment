function w = mvdr(x, theta0)
    [M, ~] = size(x);
    R = x * x' / size(x,2);

    % steering vector
    a0 = exp(1j * pi * (0:M-1)' * sind(theta0));

    % MVDR beamformer formula
    w = (R \ a0) / (a0' * (R \ a0));
end
