function [X,Y,W,SolMatT] = makeartificaldataset(N,P,K,M,nX,nY);

%make dataset
X = round(rand(N,P));


if K<=M;
    
    %randomly selecting inputs
    set = {};
    SolMatT = zeros(K,P);
    for k = 1:K;
        pos = randperm(P);
        pos = pos(1:M);
        sig = ((2*(rand(1,M)>.5)-1));
        SolMatT(k,pos) = sig;
        set{k} = pos.*sig;
    end
    
    %create binary output
    Y = zeros(N,1);
    for n = 1:length(set);
        tempvec = ones(N,1);
        for m = 1:length(set{n});
            if set{n}(m)>0;
                tempvec = tempvec&X(:,set{n}(m));
            else
                tempvec = tempvec&(~X(:,-set{n}(m)));
            end
        end
        Y = Y|tempvec;
    end
    
else
    
    %randomly selecting inputs
    set = {};
    SolMatT = zeros(M,P);
    for m = 1:M;
        pos = randperm(P);
        pos = pos(1:K);
        sig = ((2*(rand(1,K)>.5)-1));
        SolMatT(m,pos) = sig;
        set{m} = pos.*sig;
    end
    
    %create binary output
    Y = ones(N,1);
    for n = 1:length(set);
        tempvec = zeros(N,1);
        for k = 1:length(set{n});
            if set{n}(k)>0;
                tempvec = tempvec|X(:,set{n}(k));
            else
                tempvec = tempvec|(~X(:,-set{n}(k)));
            end
        end
        Y = Y&tempvec;
    end
    
end


%add noise to binary input matrix
XN = rand(size(X))<nX;
X(XN==1) = ~X(XN==1);

%add noise to binary output to create continuous output variable
Y = Y + nY*randn(size(Y));

%Create sample specific weight vector w using 0.5 as a binarization threshold
W = abs(Y-0.5);

%Binarize Y
Y = double(Y>0.5);

%Set equal class weights
W(Y==1) = W(Y==1)./sum(W(Y==1));
W(Y==0) = W(Y==0)./sum(W(Y==0));

%Set weight of samples from the Y=0 class as negative numbers
W(Y==0) = -W(Y==0);




