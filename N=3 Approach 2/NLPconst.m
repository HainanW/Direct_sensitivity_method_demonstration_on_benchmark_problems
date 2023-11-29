function [c, ceq] = NLPconst(dvar, m)

v = dvar(m.N+1:2*m.N);
c = [];
ceq = sum(v) - 8;

end 