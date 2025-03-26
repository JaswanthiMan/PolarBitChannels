clc
clear all;
close all;

N= 256;
n = log2(N);
R = 0.5;
K = floor(N*R);
mu = 64;
p_range = logspace(log10(0.05),log10(0.1),7);
BER_sim_val = zeros(1,length(p_range));
FER_sim_val = zeros(1,length(p_range));
ber_index = 1;
for p = p_range
    % erasure-method reliability sequence
    Q_em = reliability_dm(N,p,mu);
    Q_1_em = fliplr(Q_em);
    Q = Q_1_em(Q_1_em<=N);
    %Q = Q_mc(Q_mc <= N);  % erasure reliability sequence
    
    F = Q(1:N-K);
    
    Nbiterrs = 0; Nblkerrs = 0; Nblocks = 1000; Nerrs = 0;
    
    for blk = 1:Nblocks
        msg = randi([0 1],1,K); %generate random K-bit message        
        u = zeros(1,N);
        u(Q(N-K+1:end)) = msg; %assign message bits        
        cword = polar_encoding(u,n,N);    
       
        noise = (rand(1,N) <= p);%BSC channel
        y = xor(cword,noise); %received values
        
        r = get_LLR(y,N,p);
        
        msg_cap = polar_decode(r,N,F,n,K,Q);
        %Counting errors
        Nerrs = sum(msg ~= msg_cap);
        if Nerrs > 0
            Nbiterrs = Nbiterrs + Nerrs;
            Nblkerrs = Nblkerrs + 1;
        end
        %disp(blk)
    end
    BER_sim = Nbiterrs/K/Nblocks;
    FER_sim = Nblkerrs/Nblocks;

    BER_sim_val(ber_index) = BER_sim;
    FER_sim_val(ber_index) = FER_sim;
    ber_index = ber_index + 1;
    disp(p);
end



