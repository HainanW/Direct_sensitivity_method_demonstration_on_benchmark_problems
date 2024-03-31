function out = ShowMulti_Control(dvarO, m)

P = m.N;
str = cell(1,2);
str{1} = 'k-';
str{2} = 'k--';
for i = 1 : m.nu-1
    u = dvarO((i-1)*m.N + 1 : i*m.N);
    v = dvarO((m.nu-1)*m.N + 1 : end);
    ts_plot = zeros(2+2*(P-1),1);
    U_plot  = zeros(2+2*(P-1),1);
    ts_plot(1,1) = 0;
    for counter = 1 : P
        if counter == P
            ts_plot(counter*2,1) = sum(v(1:P));
             U_plot([counter*2-1 counter*2],1) = u(P);
        else
            ts_plot([2*counter 2*counter+1],1) = sum(v(1:counter));
             U_plot([counter*2-1 counter*2],1) = u(counter);
        end     
    end
    plot(ts_plot,U_plot,str{i},'LineWidth',1.5);
    hold on 
end
hold off 
xlabel('Time')
ylabel('Controls')
xlim([0 9])
% ylim([-2.83374 2.83374])
grid on 
legend({'\it{u}_{1}','\it{u}_{2}'},'Location','north')
end 