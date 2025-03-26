function [msg_cap,ucapmsg] = polar_decode(r,N,F,n,K,Q)
    % SC decoder
    L = zeros(n+1,N); %beliefs
    ucap = zeros(n+1,N); %decisions
    ns = zeros(1,2*N-1); %node state vector
    
    f = @(a,b) (1-2*(a<0)).*(1-2*(b<0)).*min(abs(a),abs(b)); %minsum
    g = @(a,b,c) b+(1-2*c).*a; %g function
    
    L(1,:) = r; %belief of root
    node = 0; depth = 0; %start at root
    done = 0; %decoder has finished or not
    while (done == 0) %traverse till all bits are decoded
        %leaf or not
        if depth == n
            if any(F==(node+1)) %is node frozen
                ucap(n+1,node+1) = 0;
            else
                if L(n+1,node+1) >= 0
                    ucap(n+1,node+1) = 0;
                else
                    ucap(n+1,node+1) = 1;
                end
            end
            if node == (N-1)
                done = 1;
            else
                node = floor(node/2); depth = depth - 1;
            end
        else
            %nonleaf
            npos = (2^depth-1) + node + 1; %position of node in node state vector
            if ns(npos) == 0 %step L and go to left child
                %disp('L')
                %disp([node depth])
                temp = 2^(n-depth);
                Ln = L(depth+1,temp*node+1:temp*(node+1)); %incoming beliefs
                a = Ln(1:temp/2); b = Ln(temp/2+1:end); %split beliefs into 2
                node = node *2; depth = depth + 1; %next node: left child
                temp = temp / 2; %incoming belief length for left child
                L(depth+1,temp*node+1:temp*(node+1)) = f(a,b); %minsum and storage
                ns(npos) = 1;
            else
                if ns(npos) == 1 %step R and go to right child
                    %disp('R')
                    %disp([node depth])
                    temp = 2^(n-depth);
                    Ln = L(depth+1,temp*node+1:temp*(node+1)); %incoming beliefs
                    a = Ln(1:temp/2); b = Ln(temp/2+1:end); %split beliefs into 2
                    lnode = 2*node; ldepth = depth + 1; %left child
                    ltemp = temp/2;
                    ucapn = ucap(ldepth+1,ltemp*lnode+1:ltemp*(lnode+1)); %incoming decisions from left child
                    node = node *2 + 1; depth = depth + 1; %next node: right child
                    temp = temp / 2; %incoming belief length for right child
                    L(depth+1,temp*node+1:temp*(node+1)) = g(a,b,ucapn); %g and storage
                    ns(npos) = 2;
                else %step U and go to parent
                    temp = 2^(n-depth);
                    lnode = 2*node; rnode = 2*node + 1; cdepth = depth + 1; %left and right child
                    ctemp = temp/2;
                    ucapl = ucap(cdepth+1,ctemp*lnode+1:ctemp*(lnode+1)); %incoming decisions from left child
                    ucapr = ucap(cdepth+1,ctemp*rnode+1:ctemp*(rnode+1)); %incoming decisions from right child
                    ucap(depth+1,temp*node+1:temp*(node+1)) = [mod(ucapl+ucapr,2) ucapr]; %combine
                    node = floor(node/2); depth = depth - 1;
                end
            end
        end
    end
    msg_cap = ucap(n+1,Q(N-K+1:end));
    ucapmsg = ucap(n+1,:);
end