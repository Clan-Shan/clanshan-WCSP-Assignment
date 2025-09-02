function [ p ] = MMSE_OSIC(Da_Str,y,H,NPW)
    F = H'*H;
    in_match = H'*y;
    p = zeros(Da_Str,1);
    dec_table = 1:Da_Str;
    for k = 1:Da_Str
        
        G = inv(F*F' + 4*NPW*F)*F;
        [~, index] = sort(diag(G'*G));
        
        w=G(:,index(1));
        z=w'*in_match;
        z = sign(real(z)) + 1i*sign(imag(z));
        t = 0.*(z == (-1+1i)) + 1.*(z == (-1-1i)) + 2.*(z == (1+1i)) + 3.*(z == (1-1i));
        
        p(dec_table(index(1))) = t;
        y=y-H(:,index(1))*(2*(real(z)>=0)-1+1i*(2*(imag(z)>=0)-1))/sqrt(2);
        
        dec_table(index(1)) = [];
        H(:,index(1)) = [];
        in_match = H'*y;
        F = H'*H;
    end
end
