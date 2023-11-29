function out = showResult_var(u,v,P)

ts_plot = zeros(2+2*(P-1),1);
U_plot  = zeros(2+2*(P-1),1);
ts_plot(1,1) = 0; 
for counter = 1: P 
    if counter == P
        ts_plot(counter*2,1) = sum(v(1:P));
        U_plot([counter*2-1 counter*2],1) = u(P);  
    else
        ts_plot([2*counter 2*counter+1],1) = sum(v(1:counter));
        U_plot([counter*2-1 counter*2],1) = u(counter);
   % tspan_real([ 2  3],1) = sum(in2(1:1));  % U_real([1 2],1) = in1(1);
   % tspan_real([ 4  5],1) = sum(in2(1:2));  % U_real([3 4],1) = in1(2);
   % tspan_real([ 6  7],1) = sum(in2(1:3));  % U_real([5 6],1) = in1(3);
   % ......
   % tspan_real( 20,1) = sum(in2(1:10)); % U_real([19 20],1) = in1(10);
        
    end
end


out = plot(ts_plot,U_plot);
xlabel('TIME')
ylabel('Control')
end 