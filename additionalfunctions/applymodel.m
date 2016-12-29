function [labels] = applymodel(x,Dm,K,M,P);

if size(Dm,2)~=P; error('Wrong number of features'); end

type = 'DNF';
if K>M;
    Ktemp = K;
    K = M;
    M = Ktemp;
    clear Ktemp;
    type = 'CNF';
end

[N,P] = size(Dm);
PM = (reshape(x(1:P*K),P,K)');
MM  = (reshape(x(P*K+1:2*P*K),P,K)');
SolMat = PM-MM;

if strcmpi(type,'DNF');
    labels = zeros(N,1);
    for k = 1:K;
        pos = find(SolMat(k,:));
        AW = SolMat(k,pos);
        AW(AW==-1) = 0;
        for m = 1:N;
            if all(AW==Dm(m,pos));
                labels(m) = 1;
            end
        end
    end
elseif strcmpi(type,'CNF');
    labels = ones(N,1);
    for k = 1:K;
        labels2 = zeros(N,1);
        pos = find(SolMat(k,:));
        AW = SolMat(k,pos);
        AW(AW==-1) = 0;
        for m = 1:N;
            if any(AW==Dm(m,pos));
                labels2(m) = 1;
            end
        end
        labels = labels&labels2;
    end
    labels = double(labels);
end

