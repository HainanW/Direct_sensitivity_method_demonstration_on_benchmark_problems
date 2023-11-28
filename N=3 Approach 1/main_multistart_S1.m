%% 
clear all 
%% Part 1.Options for ODE & NLP Solvers
odetol = 1e-8;
NLPtol = 1e-6;
optODE = odeset('RelTol', odetol, 'AbsTol', odetol);
optNLP = optimset('GradObj','on',...
    'Display', 'iter', 'TolX', NLPtol,...
    'TolFun', NLPtol, 'TolCon', NLPtol,...
    'MaxFunEvals',10000,'Algorithm','SQP');
%% Part 2.Model Pars and Algorithm Pars  
N = 3; % # of time stages 
nx = 2;
nu = 2;
nx_t = nx + nx*(nu*N); % 2+ 2*(2*3) = 14  
e1  = 18000;
e2  = 30000;
k10 = 0.535*10^11;
k20 = 0.461*10^18;
beta1 = 0.53;
beta2 = 0.43;
alpha = e2/e1;
c = k20/(k10)^alpha;

x0 = [beta1 beta2 zeros(1,nx_t-nx)]';
lb = [0.1*ones(1,N) zeros(1,N)]; 
ub = [2.0*ones(1,N) 8*ones(1,N)];
%% Run optimization with fmincon 
dvar0 = [0.166771842272956  0.144829433223544	0.129227256151918 ...
         1.63129998575985	2.49250008431194	3.87619992992821]; %JO_scaled = -0.67936690; 
dvar0 = [1*ones(1,N) 8/3*ones(1,N)];
load('functionHandle.mat')
[J0, grad0, ResVec0] = costfun1(dvar0,N,odefunvec,optODE,x0);
problem = createOptimProblem('fmincon','objective',...
    @(dvar)costfun1(dvar,N,odefunvec,optODE,x0),...
    'x0',dvar0,'options',optNLP,'lb',lb,'ub',ub,...
    'nonlcon',  @(dvar)lincon(dvar,N));
tic 
[dvarO,JO,exitflag,output,lambda,grad,hessian]  = fmincon(problem);
toc 
%% Multiple Start 
% ms = MultiStart('StartPointsToRun','bounds');
% stpoints = RandomStartPointSet('NumStartPoints',10, ...
%     'ArtificialBound',10);
% [dvarO,JO,flag,outpt,allmins] = run(ms,problem,stpoints);

%% Visualize results 
% fprintf('x2(tf)= %9.8f\n',-JO)
% u = dvarO(1:N);
% v = dvarO(N+1:2*N);
% showResult_var(u,v,N);
% % costfun1(dvarO,N,optODE)
% JO_scaled = -0.67936690; 
% dvarO =   [0.1668    0.1448    0.1292    1.6313    2.4925    3.8762];

%% Plot the discontinuity in Sensitivity 
[~, ~, ResVecO] = costfun1(dvarO,N,odefunvec,optODE,x0);
clf 
fig = figure(1);
Res = ResVecO;
t = Res(:,1);
% Sen = Res(:,10:15);
% Res: 
% Res(:,1) = t 
% 10: partial x1/t1
% 11: partial x2/t1
% 12: partial x1/t2
% 13: partial x2/t2
% 14: partial x1/t3
% 15: partial x2/t3
s11 = Res(:,10);
s21 = Res(:,11);
s12 = Res(:,12);
s22 = Res(:,13);
s13 = Res(:,14);
s23 = Res(:,15);
% plot(t, Res(:,10:15))
t_delim = [45 90];
% l1 = plt s11, s21
clf 
subplot(3,1,1)
l1 = plot(t(1:44), s11(1:44),'k--', t(1:44), s21(1:44),'k:',... 
     t(45), s11(45),'ko', t(45), s21(45),'ko', ...
     t(46), s11(46),'k.', t(46), s21(46),'k.', ...
     t(47:end), s11(47:end),'k--', t(47:end), s21(47:end), 'k:','LineWidth',1);
l1(3).MarkerSize = 4; % Make the switching point marker large as 8
l1(4).MarkerSize = 4;
l1(5).MarkerSize = 20;
l1(6).MarkerSize = 20;
xlim([0 10])
ylim([-11 7]*1e-3)
xticks(0:1:8)
% xticklabels([])
legend({'\it{\bf{\omega}}_{\it{1,1}}', ...
     '\it{\bf{\omega}}_{\it{2,1}}'},'Location','northeast',...
     'NumColumns',1);

%% 
subplot(3,1,2)
% l2 = plt s12, s22
l2 = plot(t(1:89), s12(1:89),'k--', t(1:89), s22(1:89),'k:',... 
     t(90), s12(90),'ko', t(90), s22(90),'ko', ...
     t(91), s12(91),'k.', t(91), s22(91),'k.', ...
     t(92:end), s12(92:end),'k--', t(92:end), s22(92:end), 'k:','LineWidth',1);
l2(3).MarkerSize = 4;
l2(4).MarkerSize = 4;
l2(5).MarkerSize = 20;
l2(6).MarkerSize = 20;
xlim([0 10])
ylim([-5 2.5]*1e-3)
xticks(0:1:8)
% xticklabels()
legend({'\it{\bf{\omega}}_{\it{1,2}}', ...
     '\it{\bf{\omega}}_{\it{2,2}}'},'Location','northeast',...
     'NumColumns',1);
 
subplot(3,1,3)
% l3 = plt s13, s23 
l3 = plot(t(1:134), s13(1:134),'k--', t(1:134), s23(1:134),'k:',... 
     t(end),0,'ko',t(end),0,'ko', ... 
     t(end),s13(end),'k.',t(end),s23(end),'k.','LineWidth',1);
l3(3).MarkerSize = 4;
l3(4).MarkerSize = 4;
l3(5).MarkerSize = 20;
l3(6).MarkerSize = 20;

% l4 = plot(t(45), s11(45),'k*', t(45), s21(45),'k*', ...
%      t(90), s12(90),'k*', t(90), s22(90),'k*',...
%      t(end), 0,'k*', t(end), 0,'k*');

l4(1).MarkerSize = 10;
l4(2).MarkerSize = 10;
l4(3).MarkerSize = 10;
l4(4).MarkerSize = 10;
l4(5).MarkerSize = 10;
l4(6).MarkerSize = 10;
lgd = legend({'\it{\bf{\omega}}_{\it{1,3}}', ...
    '\it{\bf{\omega}}_{\it{2,3}}'}, ... 
    'Location','northeast','NumColumns',1);
xlim([0 10])
ylim([-25 12]*1e-3)
xticks(0:1:8)
% xlabel('Time')
% ylabel('Sensitivity w.r.t. \it{t_k} ')
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Sensitivity w.r.t. \it{t_k} ');
xlabel(han,'Time (min)');


print -dsvg Figure2-L.svg