function w = mmse(x, sp)
    Rxx = x * x';        % M x M
    Rxs = x * sp';       % M x 1

    w = Rxx \ Rxs;       % Solve MMSE: w = inv(Rxx) * Rxs
end
