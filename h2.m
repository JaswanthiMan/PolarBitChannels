function h2p = h2(p)
    h2p =  -p*log2(p) - (1-p)*log2(1-p);
    if(isnan(h2p))
        h2p = 0;
    end
end