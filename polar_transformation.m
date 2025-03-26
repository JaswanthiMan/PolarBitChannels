function matrixArray = polar_transformation(Q,N)

    % Initialization of matrices with different sizes
    %N = 8;
    
    %Finding matrix sizes
    matrixSizes = zeros(N,2,log2(N)+1);
    
    for i = 1:log2(N)+1
        for j = 1:N
            matrixSizes(j, 1, i) = 2;               
            matrixSizes(j, 2, i) = 2^(2^(i-1)) * (2 ^ mod(j-1, 2^(i-1)));
        end
    end
    
    matrixArray = cell(N, log2(N)+1);
    
    for j = 1:log2(N)+1
        for i = 1:N
        matrixArray{i,j} = zeros(matrixSizes(i,1,j), matrixSizes(i,2,j));
        end
    end
    
    
    % Finding Matrix entries
    %matrixArray{1,1} = [0.2,0.8;0.8,0.2]; %Initializing first matrix
    matrixArray{1,1} = Q;
    for columns = 2:log2(N)+1
        N1 = 2^(columns-2);
        for k = 1:N1
            matrixSize = size(matrixArray{2*k-1,columns});
            W_N_i = matrixArray{k,columns-1};
            % SQUARE STAR
            for rows_k = 1:matrixSize(1)
                for columns_k = 1:matrixSize(2)
                    rows_k_bin = dec2bin(rows_k-1,log2(matrixSize(1)));
                    columns_k_bin = dec2bin(columns_k-1,log2(matrixSize(2)));
                    u2iminus1 = rows_k_bin;
                    y12N = columns_k_bin(1:2*N1);
                    y1N = y12N(1:N1);
                    yNplus12N = y12N(N1+1:end);
                    u12iminus2 = columns_k_bin(2*N1+1:end);
                    u12iminus2odd = u12iminus2(1:2:end);
                    u12iminus2even = u12iminus2(2:2:end);
                    xoru12iminus2oddeven = char((u12iminus2even ~= u12iminus2odd) + '0');
                    u2i1 = '0';
                    u2i2 = '1';
    
                    xoru2iminus1u2i1 = char((u2iminus1 ~= u2i1) + '0');
                    xoru2iminus1u2i2 = char((u2iminus1 ~= u2i2) + '0');
                    term1 = 0.5 * W_N_i(bin2dec(xoru2iminus1u2i1)+1,bin2dec([y1N,xoru12iminus2oddeven])+1)...
                            * W_N_i(bin2dec(u2i1)+1,bin2dec([yNplus12N,u12iminus2even])+1);
                    term2 = 0.5 * W_N_i(bin2dec(xoru2iminus1u2i2)+1,bin2dec([y1N,xoru12iminus2oddeven])+1)...
                            * W_N_i(bin2dec(u2i2)+1,bin2dec([yNplus12N,u12iminus2even])+1);
                    matrixArray{2*k-1,columns}(rows_k,columns_k) = term1+term2;
                end
            end
            % CIRCLE STAR
            matrixSize1 = size(matrixArray{2*k,columns});
            for rows_k1 = 1:matrixSize1(1)
                for columns_k1 = 1:matrixSize1(2)
                    rows_k1_bin = dec2bin(rows_k1-1,log2(matrixSize1(1)));
                    columns_k1_bin = dec2bin(columns_k1-1,log2(matrixSize1(2)));
                    u2i = rows_k1_bin;
                    y12N = columns_k1_bin(1:2*N1);
                    u12iminus1 = columns_k1_bin(2*N1+1:end);
                    u12iminus2 = u12iminus1(1:end-1);
                    u2iminus1 = u12iminus1(end);
                    xoru2iminus1u2i = char((u2iminus1 ~= u2i) + '0');
                    u12iminus2odd = u12iminus2(1:2:end);
                    u12iminus2even = u12iminus2(2:2:end);
                    xoru12iminus2oddeven = char((u12iminus2odd ~= u12iminus2even) + '0');
                    y1N = y12N(1:N1);
                    yNplus12N = y12N(N1+1:end);
    
                    matrixArray{2*k,columns}(rows_k1,columns_k1) = 0.5 * ...
                        W_N_i(bin2dec(xoru2iminus1u2i)+1,bin2dec([y1N,xoru12iminus2oddeven])+1) * ...
                        W_N_i(bin2dec(u2i)+1,bin2dec([yNplus12N,u12iminus2even])+1);
                end
            end
        end
    end
end
