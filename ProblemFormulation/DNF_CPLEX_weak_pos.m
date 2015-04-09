function [optfun,cons,consub,lb,ub,ctype] = DNF_CPLEX_weak_pos(X,Y,W,K,M,lambda,sens,spec,addcons)

%Quick and dirty solution to disallow negations, i.e. only positive
%predictors are allowed in the logic formulae

if nargin~=9;
    error('Please specifiy all inputs');
end

[N,P] = size(X);

%% Variables
% selector variables S and S';
% output of DNF AND gates and label for each sample (K+1)*N;
% 1 : 2*P*K + (K+1)*N

%%      Variables bound
lb = zeros(2*P*K + (K+1)*N,1);
ub = ones(2*P*K + (K+1)*N,1);

%%      Variables binary?
ctype =  char('I'*ones(1,2*P*K + (K+1)*N));

%%      Optimizing function
display('Constructing optimizing function...');
optfun = zeros(2*P*K + (K+1)*N,1);

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
cons =  sparse(NoC,2*P*K + (K+1)*N);
consub = zeros(NoC,1);
conr = 0;   %setting constraint number

%Constraint 1 OR between S and S'
display('Constructing constraint 1...');
for p = 1:P;
    for k = 1:K;
%         cons(conr+(k-1)*P+p,[(k-1)*P+p P*K+(k-1)*P+p]) = 1;
%         consub(conr+(k-1)*P+p) = 1;
        cons(conr+(k-1)*P+p,[P*K+(k-1)*P+p]) = 1;
        consub(conr+(k-1)*P+p) = 0;
    end
end
conr = P*K;

%Constraint 2 AND for kth disjuntive term
display('Constructing constraint 2...');
for n = 1:N;
    for k = 1:K;
        PosSet = find(X(n,:)==1);
        NegSet = find(X(n,:)==0);
        cons(conr+(k-1)*N+n,[2*P*K+(k-1)*N+n]) = P;
        cons(conr+(k-1)*N+n,[P*K+(k-1)*P+PosSet (k-1)*P+NegSet]) = 1;
        consub(conr+(k-1)*N+n) = P;
    end
end
conr = P*K+N*K;           

for n = 1:N;
    for k = 1:K;
        PosSet = find(X(n,:)==1);
        NegSet = find(X(n,:)==0);
        cons(conr+(k-1)*N+n,[2*P*K+(k-1)*N+n P*K+(k-1)*P+PosSet (k-1)*P+NegSet]) = -1;
        consub(conr+(k-1)*N+n) = -1;
    end
end
conr = P*K+N*K+N*K;

%Constraint 3 BIG OR between the K disjuntive terms
display('Constructing constraint 3...');
for n = 1:N;
    cons(conr+n,[2*P*K+N*K+n]) = -K;
    cons(conr+n,[2*P*K+(n:N:N*K)]) = 1;
    consub(conr+n) = 0;
end
conr = P*K+N*K+N*K+N;

for n = 1:N;
    cons(conr+n,[2*P*K+N*K+n]) = 1;
    cons(conr+n,[2*P*K+(n:N:N*K)]) = -1;
    consub(conr+n) = 0;
end
conr = P*K+N*K+N*K+N+N;


if M>0&M<P
    %Limit number of terms per disjunct
    display('Constructing constraint 4...');
    for k = 1:K;
        cons(conr+k,[(k-1)*P+1:(k-1)*P+P P*K+(k-1)*P+1:P*K+(k-1)*P+P]) = 1;
        consub(conr+k) = M;
    end
    conr = conr+K;
end

if sens>0&sens<=1;
    %Sensitivity constraint
    display('Constructing constraint 5...');
    NoP = sum(Y==1);     %Number of positives
    cons(conr+1,Io(Ip)) = -1;   %sum of TPs
    consub(conr+1) = -NoP*sens;
    conr = conr+1;
elseif sens<0&sens>=-1;
    sens = -sens;
    %Continuous sensitivity constraint
    display('Constructing constraint 5...');
    NoP = sum(W(Ip));              %Number of positives weighted by error
    cons(conr+1,Io(Ip)) = -W(Ip);  %sum of TPs
    consub(conr+1) = -NoP*sens;
    conr = conr+1;
end

if spec>0&spec<=1;
    %Specificity constraint
    display('Constructing constraint 6...');
    NoN = sum(Y==0);     %Number of negatives
    cons(conr+1,Io(In)) = 1;   %sum of FPs
    consub(conr+1) = NoN*(1-spec);
    conr = conr+1;
elseif spec<0&spec>=-1;
    spec = -spec;
    %Continuous Specificity constraint
    display('Constructing constraint 6...');
    NoN = sum(W(In));     %Number of negatives weighted by error
    cons(conr+1,Io(In)) = W(In);   %sum of FPs 
    consub(conr+1) = NoN*(1-spec);
    conr = conr+1;
end

    

if size(addcons,1)>0;
    %Additional constraints on allowed features
    display('Constructing constraint 7...');
    for c = 1:size(addcons,1);
        if strcmp(addcons{c,3},'s');
            for k = 1:K;
                cons(conr+k,[(k-1)*P+addcons{c,1} P*K+(k-1)*P+(k-1)*P+addcons{c,1}]) = 1;
                consub(conr+k) = addcons{c,2};
            end
            conr = conr+K;
        elseif strcmp(addcons{c,3},'g');
            for k = 1:K;
                cons(conr+k,[(k-1)*P+addcons{c,1} P*K+(k-1)*P+(k-1)*P+addcons{c,1}]) = -1;
                consub(conr+k) = -addcons{c,2};
            end
            conr = conr+K;
        end
    end
end




        
 