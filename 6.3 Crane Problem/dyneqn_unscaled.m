function dx = dyneqn_unscaled(t, x, u1, u2)

%DYNEQN1 
%:dynmaic state space eqn within ONE time step 

dx = zeros(7,1);
dx(1) = (x(4));
dx(2) = (x(5));
dx(3) = (x(6));
dx(4) = (u1+17.2656*x(3));
dx(5) = (u2);
dx(6) = (-(u1+27.0756*x(3)+2*x(5)*x(6))/x(2));
dx(7) = (x(3)^2+x(6)^2);

end