function [precoder] = Precoder_selection_MC(codebook, H, Da_Str, noise_power)
    codebook_size = size(codebook, 3);
    
    % Maximum Capacity Criterion
    MC_value = zeros(1, codebook_size);
    
    for index = 1:codebook_size
        % 計算當前 precoder 的 capacity
        F = codebook(:,:,index);
        capacity_matrix = eye(size(H,1)) + (1/noise_power) * H * F * (F') * (H');
        
        if det(capacity_matrix) > 0
            MC_value(1, index) = log2(det(capacity_matrix));
        else
            MC_value(1, index) = -inf; 
        end
    end
    
    % 選擇最大 capacity 對應的 precoder
    [~, F_index] = max(MC_value);
    precoder = codebook(:,:, F_index);
end