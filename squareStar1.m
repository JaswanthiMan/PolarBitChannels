function [new_Q,new_yconj] = squareStar1(Q,yconj)

    matrixSize = size(Q);
    new_matrixSize = [matrixSize(1),matrixSize(2)^2];
    new_Q = zeros(matrixSize(1),matrixSize(2)^2);

    new_yconj = zeros(1,matrixSize(2)^2);
    
    for rows = 1:new_matrixSize(1)
        rows_bin = dec2bin(rows-1,log2(new_matrixSize(1)));
        for columns = 1:new_matrixSize(2)
            columns_bin = dec2bin(columns-1,log2(new_matrixSize(2)));
            y1 = columns_bin(1:log2(matrixSize(2)));
            y2 = columns_bin(log2(matrixSize(2))+1:end);
    
            u21 = '0';
            u22 = '1';
    
            u1xoru21 = char((rows_bin ~= u21) + '0');
            u1xoru22 = char((rows_bin ~= u22) + '0');
    
            term1 = 0.5 * Q(bin2dec(u1xoru21)+1,bin2dec(y1)+1) * ...
                Q(bin2dec(u21)+1,bin2dec(y2)+1);
            term2 = 0.5 * Q(bin2dec(u1xoru22)+1,bin2dec(y1)+1) * ...
                Q(bin2dec(u22)+1,bin2dec(y2)+1);
    
            new_Q(rows,columns) = term1 + term2;

            % finding the conjugate
            columns_bin(1:log2(matrixSize(2))) = ...
                dec2bin(yconj(bin2dec(columns_bin(1:log2(...
                matrixSize(2))))+1)-1,log2(matrixSize(2)));
            new_yconj(columns) = bin2dec(columns_bin)+1;
        end
    end
end