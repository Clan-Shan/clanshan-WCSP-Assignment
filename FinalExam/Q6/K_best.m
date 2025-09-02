function [symbol_K_best] = K_best(Da_Str, y, H, Q, R)
    mod_size = 4;  % QPSK
    K = 6;         % K-best 參數
    
    sym = sqrt(1/2) * qammod(0:3, mod_size);
    
    z = Q' * y;
    
    layer = Da_Str;
    candidates = [];
    metrics = [];
    
    % 第一層：嘗試所有可能的符號
    for s = 1:mod_size
        candidate = zeros(Da_Str, 1);
        candidate(layer) = sym(s);
        
        % 計算部分距離度量
        z_est = R(layer, layer) * candidate(layer);
        metric = abs(z(layer) - z_est)^2;
        
        candidates = [candidates, candidate];
        metrics = [metrics, metric];
    end
    
    % 從第二層開始到第一層
    for layer = (Da_Str-1):-1:1
        new_candidates = [];
        new_metrics = [];
        
        % 對每個現有候選擴展
        for c = 1:size(candidates, 2)
            current_candidate = candidates(:, c);
            current_metric = metrics(c);
            
            % 嘗試當前層的所有符號
            for s = 1:mod_size
                new_candidate = current_candidate;
                new_candidate(layer) = sym(s);
                
                % 計算累積的部分距離度量
                z_est = 0;
                for j = layer:Da_Str
                    z_est = z_est + R(layer, j) * new_candidate(j);
                end
                
                partial_metric = abs(z(layer) - z_est)^2;
                total_metric = current_metric + partial_metric;
                
                new_candidates = [new_candidates, new_candidate];
                new_metrics = [new_metrics, total_metric];
            end
        end
        
        % 選擇 K 個最佳候選
        [sorted_metrics, indices] = sort(new_metrics);
        keep_num = min(K, length(indices));
        
        candidates = new_candidates(:, indices(1:keep_num));
        metrics = sorted_metrics(1:keep_num);
    end
    
    [~, best_idx] = min(metrics);
    best_symbols = candidates(:, best_idx);
    
    symbol_K_best = zeros(Da_Str, 1);
    for i = 1:Da_Str
        % 找到最接近的符號
        [~, symbol_idx] = min(abs(best_symbols(i) - sym));
        symbol_K_best(i) = symbol_idx - 1;  % 0,1,2,3
    end
end