function [new_Q, new_yconj] = circleStar1(Q,yconj)
    matrixSize = size(Q);
    new_matrixSize = [matrixSize(1),matrixSize(2)^2*2];
    new_Q = zeros(matrixSize(1),matrixSize(2)^2*2);
    
    new_yconj = zeros(1,matrixSize(2)^2*2);
    
    for rows = 1:new_matrixSize(1)
        rows_bin = dec2bin(rows-1,log2(new_matrixSize(1)));
        for columns = 1:new_matrixSize(2)
            columns_bin = dec2bin(columns-1,log2(new_matrixSize(2)));
            y1 = columns_bin(1:log2(matrixSize(2)));
            y2 = columns_bin(log2(matrixSize(2))+1:end-1);
            u1 = columns_bin(end);
    
            u1xoru2 = char((rows_bin ~= u1) + '0');
    
            new_Q(rows,columns) = 0.5*Q(bin2dec(u1xoru2)+1,bin2dec(y1)+1)*...
                Q(bin2dec(rows_bin)+1,bin2dec(y2)+1);

            % finding conjugates
            columns_bin(1:log2(matrixSize(2))) = dec2bin(yconj(bin2dec(columns_bin(1:log2(matrixSize(2))))+1)-1,log2(matrixSize(2)));
            columns_bin(log2(matrixSize(2))+1:2*log2(matrixSize(2))) = dec2bin(yconj(bin2dec(...
                columns_bin(log2(matrixSize(2))+1:2*log2(matrixSize(2))))+1)-1,log2(matrixSize(2)));
            new_yconj(columns) = bin2dec(columns_bin)+1;
    
        end
    end
end
