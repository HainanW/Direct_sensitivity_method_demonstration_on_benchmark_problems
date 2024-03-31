clear all 
clf 
close all 
%% Figure 5
load('OptRes032223_Round6.mat')
load('functionhandle.mat')
minrow_idx = 4;
dvar0_O = points(minrow_idx,:);
dvarO = dvarOMat(minrow_idx,:);
figure(1)
ShowMulti_Control(dvarO,m)
figure(2)
ShowMulti_Control(dvar0_O,m)


%% Figure 6 

tvec = dvarO(21:30);
u1vec = dvarO( 1:10);
u2vec = dvarO(11:20);

ts = zeros(1,m.N+1);
for i = 1 : m.N
    ts(i+1) = sum(tvec(1:i));
end 

z0 = x0_orig;
t_plot = [];
z_plot = [];
for i = 1 : m.N
    [t,z] = ode45(@(t,x)dyneqn_unscaled(t,x,u1vec(i),u2vec(i)),...
        [ts(i),ts(i+1)],z0,optODE);
    t_plot = [t_plot;t];
    z_plot = [z_plot;z];
    z0 = z(end,:)';
end 
figure(3)
x1_plot = z_plot(:,1);
x4_plot = z_plot(:,4);
p1 = plot(t_plot(1:10:end),x1_plot(1:10:end) ,'k-o', ... 
     t_plot,z_plot(:,2),'k:' ,...
     t_plot,z_plot(:,3),'k-',...
     t_plot(1:10:end),x4_plot(1:10:end),'k-x');
hold on 
plot(t_plot,z_plot(:,5),'k--')
plot(t_plot,z_plot(:,6),'k-.')
% plot(t_plot,z_plot(:,7),'-xk', 'MarkerSize',2)
hold off
xlim([0 9])
xlabel('Time')
ylabel('\it{x}_{\it1}~\it{x}_{6}')

legend({'\it{x}_{\it1}','\it{x}_{\it2}',...
    '\it{x}_{\it3}','\it{x}_{\it4}',...
    '\it{x}_{\it5}','\it{x}_{\it6}'},...
    'Location','Best','NumColumns',2)

%% Inequality Satisfied ? 
x4_value = z_plot(:,4);
x5_value = z_plot(:,5);

ans1 = all(abs(x4_value) -2.5)
ans2 = all(abs(x5_value) - 1)

%% Path Constraints Violation 
dvarOvec = [2.82572476918205	-0.00493465055727311	0.413989776957625	0.697068467733461	0.917876343227994	1.07388247154824	0.932626693876686	0.578397922224686	-0.109674896373806	2.83309579481293	6.36439485121201e-05	-4.24135706134010e-06	1.38694297628306e-06	2.76139178816344e-06	0.122401375399664	-0.0772732614982073	9.59148498442277e-07	0.711659418390246	0.712316099359761	0.712189821904190	0.131376530160515	0.851596182843420	1.07028685380088	1.26159274204914	1.37232918034499	2.17388501489101	0.734421222975561	0.655137298991007	0.522070148187636	0.227304825754968;
2.82572360352752	-0.00494395961209204	0.414005220318945	0.697083839429825	0.917837648975913	1.07386802380948	0.932633280167428	0.578402642826680	-0.109668107747701	2.83309555642703	5.64434225415482e-05	-5.72937874554102e-06	-1.24861646904034e-06	1.22333494758386e-05	0.122398938916631	-0.0772767738699565	-8.21039550943340e-09	0.711659212875043	0.712315898059281	0.712189599532673	0.131391197295863	0.851608650618369	1.07028664443622	1.26158796608500	1.37232537575052	2.17384153527360	0.734443321142236	0.655124768397732	0.522069347772786	0.227321193227700];
sigma = 1e-6;
rho = 1;
tf = 10;
for i = 1 : 2
    rho = 10^(i-1);
    sigma = 1e-6;
    epsilon = 2*sigma/((1+sqrt(5))*rho*2*tf);
    costfunTS(dvarOvec(i,:), m, odefunvec, optODE, x0_red, epsilon, rho)
end 