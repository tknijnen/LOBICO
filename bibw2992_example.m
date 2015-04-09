%% Initialize
Init

%% Load data

%parse tsv
% [Samples,Features,IC50s,MutationMatrix] = parsetsv('bibw2992.tsv');

%parse xls file
% [Samples,Features,IC50s,MutationMatrix] = parsexls('bibw2992.xls');

%load matfile
load bibw2992

%% Create binary input, output and weight vector

%binary input
X = MutationMatrix;
[N,P] = size(X);

%binarization threshold th
th = 0.063218; 
Y = double(IC50s<th);
W = abs(IC50s-th);

%class weights
FPW = 1;                  %Total weight on positive class (Y==1)
FPN = 1;                  %Total weight on negative class (Y==0)

%normalize weights
W(Y==1) = FPW*W(Y==1)./sum(W(Y==1));
W(Y~=1) = -(FPN*W(Y~=1)./sum(W(Y~=1)));

%% Logic model complexity
K = 2;
M = 1;

%% Cplex options
param = cat(1,{'timelimit.Cur',60,'MaxTime'},...                            %Maximum time for IP (in seconds)
              {'mip.tolerances.integrality.Cur',1e-5,'Integrality'},...     %Integrality contraint; default = 1e-5 (see cplex.Param.mip.tolerances.integrality.Help)
              {'mip.tolerances.mipgap.Cur',1e-4,'RelGap'},...               %Optimality gap tolerance; default = 1e-4 (0.01% of optimal solution, set to 0.05, 5% for early termination, approximate solution) (see cplex.Param.mip.tolerances.mipgap.Help)
              {'threads.Cur',8,'Threads'},...                               %Number of threads to use (default = 0, automatic) (see  cplex.Param.threads.Help);          
              {'parallel.Cur',-1,'ParallelMode'},...                        %Parallel optimization mode,  -1 = opportunistic 0 = automatic 1 = deterministic (see cplex.Param.parallel.Help)
              {'mip.pool.relgap.Cur',1e-1,'Pool_relgap'},...                %Relative gap for suboptimal solutions in the pool; default 0.1 (10%)
              {'mip.pool.intensity.Cur',1,'Pool_intensity'},...             %Pool intensity; default 1 = mild: generate few solutions quickly (see  cplex.Param.mip.pool.intensity.Help)
              {'mip.pool.replace.Cur',2,'Pool_replace'},...                 %Pool replacement strategy; default 2 = replace least diverse solutions (see  cplex.Param.mip.pool.replace.Help)
              {'mip.pool.capacity.Cur',11,'Pool_capacity'},...              %Pool capacity; default 11 = best + 10 (see  cplex.Param.mip.pool.replace.Help)
              {'mip.limits.populate.Cur',11,'Pool_size'});                  %Number of solutions generated; default 11 = best + 10 (see  cplex.Param.mip.limits.populate.Help)


%% Cplex solver
sol = lobico(X,W,K,M,1,param);

%% Check solution
display('***********************');

%inferred formula
x = round(sol.Solution.x);
SolMat = getsolution(x,K,M,P);
str = showformula(SolMat,K,M,Features);
disp('Inferred logic model');
disp(str);

display('***********************');
