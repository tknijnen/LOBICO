function str = showformula(SolMat,K,M,Xl)

if ~exist('Xl','var');
    Xl = cell(size(SolMat,2),1);
    for n = 1:length(Xl);
        Xl{n} = ['x' num2str(n)];
    end
end
    

if K<=M;

    str = [];
    for k = 1:K;
        pos = find(SolMat(k,:));
        for p = 1:length(pos);
            if SolMat(k,pos(p))==-1;
                str = [str '~'];
            end
            str = [str Xl{pos(p)} ' '];
            if p~=length(pos);
                str = [str ' & '];
            else
                if k~=K;
                    str = [str '  |  '];
                end
            end
        end
    end
    
else
  
    Ktemp = K;
    K = M;
    M = Ktemp;
    
    str = [];
    for k = 1:K;
        pos = find(SolMat(k,:));
        for p = 1:length(pos);
            if SolMat(k,pos(p))==-1;
                str = [str '~'];
            end
            str = [str Xl{pos(p)} ' '];
            if p~=length(pos);
                str = [str ' | '];
            else
                if k~=K;
                    str = [str '  &  '];
                end
            end
        end
    end
end