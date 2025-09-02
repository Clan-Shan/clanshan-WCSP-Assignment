function [decide_data] = MMSE_OSIC_2025(H, in, P, NPW, M, L)
% H : is the channel matrix
% in : is the input data
% P : is the number of data streams
% NPW : is the noise power
% 4QAM used: M=4, L=2
F = H'*H; % Gram 矩陣，用於計算 MMSE 濾波器
in_match = H'*in;
decide_data = zeros(1,P); % 存放解碼結果
dec_table = 1:P; % a index table helps make decision in tx_order

for k = 1:P
    G = inv(F*F' + NPW*F)*F;
    [~, index] = sort(diag(G'*G));
    
    Z=G(:,index(1))'*in_match;
    decide = qamdemod(Z,M);  % QAM 解調
    decide_data(dec_table(index(1))) = decide;
    dec = qammod(decide,M); % 再調變回去，用來重建訊號
    in=in-H(:,index(1))*1/(L^0.5)*dec; % 消除此符號對整體的貢獻
    
    dec_table(index(1)) = [];%移除該欄，繼續解剩下的符號
    H(:,index(1)) = [];
    in_match = H'*in;
    F = H'*H;
end

end