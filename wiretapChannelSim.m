clc
%wire-tap construction.

clear all;
close all;

%Eve Information
EveMatrix = readmatrix('IPZ_Eve_015_256.txt');
%Bob Information
BobMatrix = readmatrix('IPZ_Bob_005.txt');

IEve = EveMatrix(1,:);
PEve = EveMatrix(2,:);
ZEve = EveMatrix(3,:);

IBob = BobMatrix(1,:);
PBob = BobMatrix(2,:);
ZBob = BobMatrix(3,:);

N = 256;
n = log2(N);
delta = 0.25;
beta = 0.15;
beta1 = 0.2;
beta_n = 2^(-N^(beta));
beta_n1 = 2^(-N^(beta1));
p2 = 0.4;
p1 = 0.05;
capacity = h2(p2) - h2(p1);


[P_Bob,QBob] = sort(PBob);
%[Z_sorted,Q_em] = reliability_em(N,p);
Q = fliplr(QBob);
ppp = find((fliplr(P_Bob) <= 2^(-N^(beta))));
Gbob = Q(ppp); % Good set of Bob
GBobcomp = setdiff(Q,Gbob); 

%%
[C_Eve,QEve] = sort(IEve);
Qeve = fliplr(QEve);
pppeve = find((fliplr(C_Eve) < 0.3));
PoorEve = Qeve(pppeve);
PoorEvecomp = setdiff(Qeve,PoorEve);
setA = intersect(Gbob,PoorEve);
setB = intersect(PoorEve,GBobcomp);
setX = intersect(GBobcomp,PoorEvecomp);
SetY = intersect(Gbob,PoorEvecomp);

keycapacity = h2(p2) - h2(p1);

%%

key_Rate = length(setA)/N;

%%
keylength = length(setA);
F = setB;
R = setX;
K = N - length(setB);
Nblocks = 10000; Nblockerrsbob = 0; Nblockerrseve=0;
for i = 1:Nblocks
    V = randi([0,1],1,N);
    V(setB) = 0; 
    cword = polar_encoding(V,n,N);
    %channel
    noise = (rand(1,N) <= p1);%BSC channel
    Y = xor(cword,noise); %received value
    %receiver
    r = get_LLR(Y,N,p1);
    [msg_cap,Vcap] = polar_decode_Bob(r,N,F,n,K,Q,R,V);
    % Eve Decoding
    noise = (rand(1,N) <= p2);
    Z = xor(cword,noise);
    reve = get_LLR(Z,N,p2);
    [msg_capeve,Vcapeve] = polar_decode(reve,N,F,n,K,Q);
    if(~isequal(Vcap(setA), V(setA)))
        Nblockerrsbob = Nblockerrsbob + 1;
    end
    if(~isequal(Vcapeve(setA), V(setA)))
        Nblockerrseve = Nblockerrseve + 1;
    end
    %Eve Maximum likelihood decoding
    Receivedbits = Z(setdiff(1:N,setB));


end
% 
% %%
% % keyrate = [0.3008,0.2578,0.2305,0.1289,0.0742,0.0430];
% % BobError = [0.4727,0.3688,0.1782,0.0231,0.0083,0.0062];
% % beta = [0.1,0.12,0.15,0.17,0.19,0.22,0.26];
% % 
% % %Bob FER
% % loglog(keyrate,BobError,'-*');
% % xlabel('Key Rate');
% % ylabel('Bob FER');
% % title('Bob FER at key capacity = 0.5949,N = 256,p = 0.05,q=0.3');
% 
% % At length N = 64
% 
% % keyrate = [0.1875,0.1094,0.0625,0.0313,0.0156];
% % beta = [0.18,0.25,0.3,0.33,0.38];
% % Bob_FER = [0.2970,0.1348,0.0211,0.0111,0.0064];
% % Eve_FER = [0.9997,0.9897,0.9297,0.7458,0.4922];
% % 
% % loglog(keyrate,Bob_FER,'-*');
% % xlabel('Key Rate');
% % ylabel('Bob FER');
% % title('Bob FER at key capacity = 0.5949,N = 64,p = 0.05,q=0.3,mu=32');
% 
% %%
% % Implementation of maximum likelihood decoder for Eve
% 

%%

% clear all;
% p = 0.05;
% q = [0.15,0.2,0.25,0.3,0.35,0.4];
% keycapacity = zeros(1,length(q));
% for i = 1:length(q)
%     keycapacity(i) = h2(q(i))-h2(p);
% end
% achievablerate = [0.29,0.34,0.4,0.44,0.46,0.48];
% actualrate = [0.08,0.2,0.3,0.35,0.4,0.45];
% 
% figure;
% plot(q,keycapacity, '*-');
% hold on;
% plot(q, achievablerate, '*-');
% hold on;
% plot(q, actualrate, '*-');
% xlabel('q');
% ylabel('Rate');
% legend('Secure Message Cap', 'Achievable Rate SCD', 'Actual Secure Message Rate');
% title('Secure Message Protocol at p = 0.05, N = 256 while maintaining Bob FER~0.5');