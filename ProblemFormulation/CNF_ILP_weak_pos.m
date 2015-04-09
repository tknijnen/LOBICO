function [optfun,cons,consub,lb,ub,ctype] = CNF_ILP_weak_pos(X,Y,W,K,M,lambda,sens,spec,addcons)

%Quick and dirty solution to disallow negations, i.e. only positive
%predictors are allowed in the logic formulae

if nargin~=9;
    error('Please specifiy all inputs');
end

[N,P] = size(X);

%% Variables
% selector variables S and S';
% output of CNF OR gates and label for each sample (K+1)*N;
NoV = 2*P*K + (K+1)*N;

%%      Variables bound
lb = zeros(NoV,1);
ub = ones(NoV,1);

%%      Variables binary?
ctype =  char('I'*ones(1,NoV));

%%      Optimizing function
display('Constructing optimizing function...');
optfun = zeros(NoV,1);

%error
Ip = find(Y==1);
In = find(Y==0);
Io = 2*P*K + K*N + 1:2*P*K + (K+1)*N;
optfun(Io(In)) = W(In);                                 %Error if mispredicting a 0
optfun(Io(Ip)) = -W(Ip);                                %Error if mispredicting a 1

%model complexity penalty
if M==-1;
    optfun(1:2*P*K) = lambda;
end

%%      Constraints and their bounds
display('Constructing constraints...');
%number of constraints
NoC = P*K+N*K+N*K+N+N;
if M>0&M<P
    NoC = NoC+K;
end
if (sens>0&sens<=1)|(sens<0&sens>=-1);
    NoC = NoC+1;
end
if (spec>0&spec<=1)|(spec<0&spec>=-1);
    NoC = NoC+1;
end
NoC = NoC + K*size(addcons,1);

%number of nonzero elements
NoZ = P*K+N*K+K*N*P+N*K+K*N*P+N+N*K+N+N*K;
if M>0&M<P
    NoZ = NoZ+2*P*K;
end
if (sens>0&sens<=1)|(sens<0&sens>=-1);
    NoZ = NoZ+length(Ip);
end
if (spec>0&spec<=1)|(spec<0&spec>=-1);
    NoZ = NoZ+length(In);
end
for c = 1:size(addcons,1);
    NoZ = NoZ+K*2*length(addcons{c,1});
end

%make index variables for sparse matrix
I = zeros(NoZ,1);
J = zeros(NoZ,1);
Q = zeros(NoZ,1);

consub = zeros(NoC,1);
conr = 0;   %setting constraint number
conz = 0;   %setting number of nonzero elements

%Constraint 1 OR between S and S'
p = 1:P;
for k = 1:K;
    ix = conz+p;
    I(ix) = conr+(k-1)*P+p;
    J(ix) = P*K+(k-1)*P+p;
    Q(ix) = 1;
    consub(conr+(k-1)*P+p) = 0;
    conz = conz+P;
end
conr = P*K;
conz = P*K;

%Constraint 2 OR for kth conjuntive term
k = 1:K;
for n = 1:N;
    PosSet = find(X(n,:)==1);
    NegSet = find(X(n,:)==0);
    LPS = length(PosSet);
    LNS = length(NegSet);
    ix = conz+k;
    I(ix) = conr+(k-1)*N+n;
    J(ix) = 2*P*K+(k-1)*N+n;
    Q(ix) = 1;
    conz = conz+K;
    ix = conz+k;
    I(ix) = conr+N*K+(k-1)*N+n;
    J(ix) = 2*P*K+(k-1)*N+n;
    Q(ix) = -P;
    conz = conz+K;
    for kk=1:K;
        ix = conz+(1:LPS);
        I(ix) = conr+(kk-1)*N+n;
        J(ix) = (kk-1)*P+PosSet;
        Q(ix) = -1;
        conz = conz+LPS;
        ix = conz+(1:LNS);
        I(ix) = conr+(kk-1)*N+n;
        J(ix) = P*K+(kk-1)*P+NegSet;
        Q(ix) = -1;
        conz = conz+LNS;
        ix = conz+(1:LPS);
        I(ix) = conr+N*K+(kk-1)*N+n;
        J(ix) = (kk-1)*P+PosSet;
        Q(ix) = 1;
        conz = conz+LPS;
        ix = conz+(1:LNS);
        I(ix) = conr+N*K+(kk-1)*N+n;
        J(ix) = P*K+(kk-1)*P+NegSet;
        Q(ix) = 1;
        conz = conz+LNS;
    end
    consub(conr+(k-1)*N+n) = 0;
    consub(conr+N*K+(k-1)*N+n) = 0;
end
conr = P*K+N*K+N*K;
conz = P*K+N*K+K*N*P+N*K+K*N*P;

%Constraint 3 BIG AND between the K conjuntive terms
n = 1:N;
ix = conz+n;
I(ix) = conr+n;
J(ix) = 2*P*K+N*K+n;
Q(ix) = -1;
conz = conz+N;
ix = conz+n;
I(ix) = conr+N+n;
J(ix) = 2*P*K+N*K+n;
Q(ix) = K;
conz = conz+N;
k=1:K;
for nn=1:N;
    ix = conz+k;
    I(ix) = conr+nn;
    J(ix) = 2*P*K+(nn:N:N*K);
    Q(ix) = 1;
    conz = conz+K;
    ix = conz+k;
    I(ix) = conr+N+nn;
    J(ix) = 2*P*K+(nn:N:N*K);
    Q(ix) = -1;
    conz = conz+K;
end
consub(conr+n) = K-1;
consub(conr+N+n) = 0;
conr = P*K+N*K+N*K+N+N;
conz = P*K+N*K+K*N*P+N*K+K*N*P+N+N*K+N+N*K;

if M>0&M<P
    %Limit number of terms per disjunct
    p = 1:P;
    for k = 1:K;
        ix = conz+p;
        I(ix) = conr+k;
        J(ix) = (k-1)*P+1:(k-1)*P+P;
        Q(ix) = 1;
        conz = conz+P;
        ix = conz+p;
        I(ix) = conr+k;
        J(ix) = P*K+(k-1)*P+1:P*K+(k-1)*P+P;
        Q(ix) = 1;
        conz = conz+P;
        consub(conr+k) = M;
    end
    conr = conr+K;
end

if sens>0&sens<=1;
    %Sensitivity constraint
    NoP = sum(Y==1);     %Number of positives
    ix = conz+(1:length(Ip));
    I(ix) = conr+1;
    J(ix) = Io(Ip);
    Q(ix) = -1;
    conz = conz+length(Ip);
    consub(conr+1) = -NoP*sens;
    conr = conr+1;
elseif sens<0&sens>=-1;
    sens = -sens;
    %Continuous sensitivity constraint
    NoP = sum(W(Ip));              %Number of positives weighted by error
    ix = conz+(1:length(Ip));
    I(ix) = conr+1;
    J(ix) = Io(Ip);
    Q(ix) = -W(Ip);
    conz = conz+length(Ip);
    consub(conr+1) = -NoP*sens;
    conr = conr+1;
end

if spec>0&spec<=1;
    %Specificity constraint
    NoN = sum(Y==0);     %Number of negatives
    ix = conz+(1:length(In));
    I(ix) = conr+1;
    J(ix) = Io(In);
    Q(ix) = 1;
    conz = conz+length(In);
    consub(conr+1) = NoN*(1-spec);
    conr = conr+1;
elseif spec<0&spec>=-1;
    spec = -spec;
    %Continuous Specificity constraint
    NoN = sum(W(In));     %Number of negatives weighted by error
    ix = conz+(1:length(In));
    I(ix) = conr+1;
    J(ix) = Io(In);
    Q(ix) = W(In);
    conz = conz+length(In);
    consub(conr+1) = NoN*(1-spec);
    conr = conr+1;
end

if size(addcons,1)>0;
    %Additional constraints on allowed features
    for c = 1:size(addcons,1);
        if strcmp(addcons{c,3},'s');
            for k = 1:K;
                ix = conz+(1:length(addcons{c,1}));
                I(ix) = conr+k;
                J(ix) = (k-1)*P+addcons{c,1};
                Q(ix) = 1;
                conz = conz+length(addcons{c,1});
                ix = conz+(1:length(addcons{c,1}));
                I(ix) = conr+k;
                J(ix) = P*K+(k-1)*P+(k-1)*P+addcons{c,1};
                Q(ix) = 1;
                conz = conz+length(addcons{c,1});
                consub(conr+k) = addcons{c,2};
            end
            conr = conr+K;
        elseif strcmp(addcons{c,3},'g');
            for k = 1:K;
                ix = conz+(1:length(addcons{c,1}));
                I(ix) = conr+k;
                J(ix) = (k-1)*P+addcons{c,1};
                Q(ix) = -1;
                conz = conz+length(addcons{c,1});
                ix = conz+(1:length(addcons{c,1}));
                I(ix) = conr+k;
                J(ix) = P*K+(k-1)*P+(k-1)*P+addcons{c,1};
                Q(ix) = -1;
                conz = conz+length(addcons{c,1});
                consub(conr+k) = -addcons{c,2};
            end
            conr = conr+K;
        end
    end
end

cons =  sparse(I,J,Q,NoC,NoV,NoZ);
    
    
    
    
    
 