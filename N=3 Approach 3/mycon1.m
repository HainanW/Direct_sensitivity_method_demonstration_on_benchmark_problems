function [c, ceq] = mycon1(dvar, N)

v = dvar(N+1:2*N);
c = [];
ceq = sum(v) - 8;

end 