function mutual_information_W = compute_mutualinformation(W)
    size_W = size(W);
    num_columns = size_W(2);
    num_rows = size_W(1);
    mutual_information_W = 0;
    for y = 1:num_columns
        row_I = 0;
        for x = 1:num_rows
            row_I = row_I + 0.5*W(1,y)*log2(W(1,y)/(0.5*W(1,y)+0.5*W(2,y)));
        end
        mutual_information_W = mutual_information_W + row_I;
    end
end