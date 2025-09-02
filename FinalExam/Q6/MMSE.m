function [ p ] = MMSE(Da_Str,y,H,NPW)
    p = zeros(Da_Str,1);
    x_hat = inv(H'*H + NPW * eye(Da_Str)) * H' * y;

    p = qamdemod(x_hat,4);
end
