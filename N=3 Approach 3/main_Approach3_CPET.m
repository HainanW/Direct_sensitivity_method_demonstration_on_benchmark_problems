clf
clear all 
%% Part 1.Options for ODE & NLP Solvers
odetol = 1e-8;
NLPtol = 1e-6;
optODE = odeset('RelTol', odetol, 'AbsTol', odetol);
optNLP = optimset('GradObj','on', 'GradConstr','off',...
    'DerivativeCheck', 'off', 'Display', 'iter', 'TolX', NLPtol,...
    'TolFun', NLPtol, 'TolCon', NLPtol, 'MaxFunEvals',10000,'Algorithm','SQP');
%% Part 2.Model Pars and Algorithm Pars  
N = 3; % # of time stages 
nx = 2;
nu = 2;
nx_t = nx + nx*(nu*N); % 2+ 2*(2*6) = 26 
e1  = 18000;
e2  = 30000;
k10 = 0.535*10^11;
k20 = 0.461*10^18;
beta1 = 0.53;
beta2 = 0.43;
alpha = e2/e1;
c = k20/(k10)^alpha;
tf = N;
ts = linspace(0,tf,N+1);
x0 = [beta1 beta2 zeros(1,nx_t-nx)]';
lb = [0.1*ones(1,N) zeros(1,N)]; 
ub = [2.0*ones(1,N) Inf*ones(1,N)];
%% Part 3. Run optimization with fmincon 
load('functionHandle.mat')
dvar0 = [0.166771842272956	0.144829433223544	0.129227256151918 ...
         1.63129998575985	2.49250008431194	3.87619992992821]; %JO_scaled = -0.67936690; 
[J0, grad0, zvec_scaled] = costfun1(dvar0,N,odefunvec,optODE,ts,x0);
dvar0 = [1*ones(1,N) 8/3*ones(1,N)];
problem = createOptimProblem('fmincon','objective',...
    @(dvar)costfun1(dvar,N,odefunvec,optODE,ts,x0),...
    'x0',dvar0,'options',optNLP,'lb',lb,'ub',ub,...
    'nonlcon',  @(dvar)mycon1(dvar,N));
tic 
[dvarO,JO] = fmincon(problem);
toc


%% Visualize results 

fprintf('x2(tf)= %9.8f\n',-JO)
[JO, gradO, zvec_scaled] = costfun1(dvarO,N,odefunvec,optODE,ts,x0);
u = dvarO(1:N);
v = dvarO(N+1:2*N);
% showResult_var(u,v,N);
save('zvec_scaled.mat','zvec_scaled')