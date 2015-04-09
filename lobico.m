function cplex = lobico(X,Y,K,M,solve,param,spec,sens,lambda,weak,pos,addcons);

%% Set undefined input arguments to default settings
if ~exist('solve','var');    solve   = 1;           elseif isempty(solve);     solve   = 1;         end
if ~exist('param','var');    param  = cell(0,3);    elseif isempty(param);     param  = cell(0,3);  end
if ~exist('spec','var');     spec   = 0;            elseif isempty(spec);      spec   = 0;          end
if ~exist('sens','var');     sens   = 0;            elseif isempty(sens);      sens   = 0;          end
if ~exist('lambda','var');   lambda = 0;            elseif isempty(lambda);    lambda = 0;          end
if ~exist('weak','var');     weak   = 1;            elseif isempty(weak);      weak   = 1;          end
if ~exist('pos','var');      pos    = 0;            elseif isempty(pos);       pos    = 0;          end
if ~exist('addcons','var');  addcons= cell(0,2);    elseif isempty(addcons);   addcons= cell(0,2);  end

%% Formulate problem
if K<=M;
    cplex = Cplex('DNF');
    if pos==1
        [f,Aineq,bineq,lb,ub,ctype] = DNF_ILP_weak_pos(X,double(Y>0),abs(Y),K,M,lambda,sens,spec,addcons);
    else
        if weak==1
            [f,Aineq,bineq,lb,ub,ctype] = DNF_ILP_weak(X,double(Y>0),abs(Y),K,M,lambda,sens,spec,addcons);
        else
            [f,Aineq,bineq,lb,ub,ctype] = DNF_CPLEX(X,double(Y>0),abs(Y),K,M,lambda,sens,spec,addcons);
        end
    end
else
    cplex = Cplex('CNF');
    if pos==1
        [f,Aineq,bineq,lb,ub,ctype] = CNF_ILP_weak_pos(X,double(Y>0),abs(Y),M,K,lambda,sens,spec,addcons);
    else
        if weak==1;
            [f,Aineq,bineq,lb,ub,ctype] = CNF_ILP_weak(X,double(Y>0),abs(Y),M,K,lambda,sens,spec,addcons);
        else
            [f,Aineq,bineq,lb,ub,ctype] = CNF_CPLEX(X,double(Y>0),abs(Y),M,K,lambda,sens,spec,addcons);
        end
    end
end
cplex.Model.obj   = f;
cplex.Model.sense = 'minimize';
cplex.Model.lb = lb;
cplex.Model.ub = ub;
cplex.Model.ctype = ctype;
cplex.Model.A   = Aineq;
cplex.Model.rhs = bineq;
cplex.Model.lhs = -Inf*ones(size(bineq));

%% Set parameters
try
    np = size(param,1);
    for n = 1:np;
        try
            eval(['cplex.Param.' param{n,1} ' = param{n,2};']);
            disp([param{n,3} ' set to ' num2str(param{n,2})]);
        catch
            disp(lasterr);
            disp(['Setting ' param{n,3} ' failed']);
        end
    end
catch
    disp(['Setting parameters failed'])
end

%% Solving problem
if solve==1;
    cplex.solve;
elseif solve==2
    cplex.populate;
end
