function SolMat = getsolution(x,K,M,P);

if K>M;
    Ktemp = K;
    K = M;
    M = Ktemp;
    clear Ktemp;
end

PM = (reshape(x(1:P*K),P,K)');
MM  = (reshape(x(P*K+1:2*P*K),P,K)');
SolMat = PM-MM;