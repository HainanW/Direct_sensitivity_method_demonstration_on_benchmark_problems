%%% Trajectory Sensitivity Approach; 
%%% Variable Length 
%% 052621 Result Record 
dvarO = [0.174175775948736,0.161006422056221,0.150014493556339,0.140768853991075,0.132953496404172,0.126334773741388,0.728008446146651,1.62939081472582,2.74296366962682,4.12098873365750,5.83695915540293];
u = dvarO(1:6);
z0 = [0.53,0.43]';
ts = linspace(0,8,7);
ts(7) = 8;
ts(2:6) = dvarO(7:11);
for ks = 1 : 6
    [~,zs] = ode15s(@(t,x)dyneqn(t,x,u(ks)), ...
        [ts(ks),ts(ks+1)], z0, optODE);
    z0 = zs(end,:)';
end
f = zs(end,:)';
