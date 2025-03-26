function r = get_LLR(y,N,p)
    r = zeros(1,N); %beliefs
    %finding the belief vector corresponding to the received codeword
    for i = 1:N
        if y(i) == 0
            r(i) = log((1-p)/p);
        else
            r(i) = log(p/(1-p));
        end
    end
end