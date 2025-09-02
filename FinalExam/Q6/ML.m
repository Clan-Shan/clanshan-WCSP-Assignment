function [symbol_ML] = ML(Da_Str, y, H, Q, R, NPW)
    % 窮舉搜索所有可能的符號組合
    
    symbol_ML = zeros(Da_Str, 1);
    
    symbol_d = (1/sqrt(2)) * [1+1j, 1-1j, -1+1j, -1-1j];
    
    if Da_Str == 1
        all_possible_x = symbol_d.';
    elseif Da_Str == 2
        all_possible_x = combvec(symbol_d, symbol_d);
    elseif Da_Str == 3
        all_possible_x = combvec(symbol_d, symbol_d, symbol_d);
    elseif Da_Str == 4
        all_possible_x = combvec(symbol_d, symbol_d, symbol_d, symbol_d);
    else
        all_possible_x = generate_all_combinations(symbol_d, Da_Str);
    end
    
    minDist = inf;
    detectedSymbolIdx = 1;
    [~, m] = size(all_possible_x);
    
    for symbol = 1:m
        candidate = all_possible_x(:, symbol);
        
        % 計算距離
        dist = norm(y - H * candidate)^2;
        
        if dist < minDist
            minDist = dist;
            detectedSymbolIdx = symbol;
        end
    end
    
    detected_symbols = all_possible_x(:, detectedSymbolIdx);
    
    for i = 1:Da_Str
        symbol_ML(i) = qpsk_symbol_to_index(detected_symbols(i));
    end
end

function all_combinations = generate_all_combinations(symbol_set, num_streams)
    if num_streams == 1
        all_combinations = symbol_set.';
    else
        prev_combinations = generate_all_combinations(symbol_set, num_streams - 1);
        num_symbols = length(symbol_set);
        num_prev = size(prev_combinations, 2);
        
        all_combinations = zeros(num_streams, num_symbols * num_prev);
        idx = 1;
        
        for s = 1:num_symbols
            for p = 1:num_prev
                all_combinations(:, idx) = [symbol_set(s); prev_combinations(:, p)];
                idx = idx + 1;
            end
        end
    end
end

function idx = qpsk_symbol_to_index(symbol)
    %  0→-1+1j, 1→-1-1j, 2→1+1j, 3→1-1j
    
    symbol_scaled = symbol * sqrt(2);  
    
    if abs(symbol_scaled - (1+1j)) < 1e-10
        idx = 2;  % 1+1j
    elseif abs(symbol_scaled - (1-1j)) < 1e-10
        idx = 3;  % 1-1j
    elseif abs(symbol_scaled - (-1+1j)) < 1e-10
        idx = 0;  % -1+1j
    elseif abs(symbol_scaled - (-1-1j)) < 1e-10
        idx = 1;  % -1-1j
    else
        % 找最接近的符號
        qpsk_symbols = [(-1+1j), (-1-1j), (1+1j), (1-1j)];
        [~, idx] = min(abs(symbol_scaled - qpsk_symbols));
        idx = idx - 1;  
    end
end