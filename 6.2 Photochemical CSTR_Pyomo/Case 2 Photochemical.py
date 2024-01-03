import idaes
from pyomo.environ import *
from pyomo.dae import *
from idaes.core.util.model_diagnostics import DegeneracyHunter

# Create a model project
m = ConcreteModel()
N = 66
m.n = RangeSet(N)
tauF = 0.2 / N
m.tau = ContinuousSet(bounds=(0, tauF))
# Define the Control variables

m.u1 = Var(m.n, bounds=(0, 20), initialize=10)
m.u2 = Var(m.n, bounds=(0, 6), initialize=3)
m.u3 = Var(m.n, bounds=(0, 4), initialize=2)
m.u4 = Var(m.n, bounds=(0, 20), initialize=10)
m.q = Var(m.n, domain=NonNegativeReals)
m.qexpr = Constraint(m.n, rule=lambda m, n: m.q[n] == m.u1[n] + m.u2[n] + m.u4[n])

m.x1 = Var(m.n, m.tau, domain=NonNegativeReals, initialize=0.1883)
m.x2 = Var(m.n, m.tau, domain=NonNegativeReals, initialize=0.2507)
m.x3 = Var(m.n, m.tau, domain=NonNegativeReals, initialize=0.0467)
m.x4 = Var(m.n, m.tau, domain=NonNegativeReals, initialize=0.0899)
m.x5 = Var(m.n, m.tau, domain=NonNegativeReals, initialize=0.1804)
m.x6 = Var(m.n, m.tau, domain=NonNegativeReals, initialize=0.1394)
m.x7 = Var(m.n, m.tau, domain=NonNegativeReals, initialize=0.1046)
m.x8 = Var(m.n, m.tau, domain=NonNegativeReals, initialize=0.0000)

m.dx1 = DerivativeVar(m.x1, wrt=m.tau)
m.dx2 = DerivativeVar(m.x2, wrt=m.tau)
m.dx3 = DerivativeVar(m.x3, wrt=m.tau)
m.dx4 = DerivativeVar(m.x4, wrt=m.tau)
m.dx5 = DerivativeVar(m.x5, wrt=m.tau)
m.dx6 = DerivativeVar(m.x6, wrt=m.tau)
m.dx7 = DerivativeVar(m.x7, wrt=m.tau)
m.dx8 = DerivativeVar(m.x8, wrt=m.tau)

m.ode1 = Constraint(m.n, m.tau,
                    rule=lambda m, n, t: m.dx1[n, t] == m.u4[n] - m.q[n] * m.x1[n, t] - 17.6 * m.x1[n, t] * m.x2[n, t] -
                                         23.0 * m.x1[n, t] * m.x6[n, t] * m.u3[n])

m.ode2 = Constraint(m.n, m.tau,
                    rule=lambda m, n, t: m.dx2[n, t] == m.u1[n] - m.q[n] * m.x2[n, t] - 17.6 * m.x1[n, t] * m.x2[n, t] -
                                         146 * m.x2[n, t] * m.x3[n, t])

m.ode3 = Constraint(m.n, m.tau,
                    rule=lambda m, n, t: m.dx3[n, t] == m.u2[n] - m.q[n] * m.x3[n, t] - 73 * m.x2[n, t] * m.x3[n, t])

m.ode4 = Constraint(m.n, m.tau,
                    rule=lambda m, n, t: m.dx4[n, t] == - m.q[n] * m.x4[n, t] + 35.2 * m.x1[n, t] * m.x2[n, t]
                                         - 51.3 * m.x4[n, t] * m.x5[n, t])

m.ode5 = Constraint(m.n, m.tau,
                    rule=lambda m, n, t: m.dx5[n, t] == - m.q[n] * m.x5[n, t] + 219 * m.x2[n, t] * m.x3[n, t]
                                         - 51.3 * m.x4[n, t] * m.x5[n, t])

m.ode6 = Constraint(m.n, m.tau,
                    rule=lambda m, n, t: m.dx6[n, t] == - m.q[n] * m.x6[n, t] + 102.6 * m.x4[n, t] * m.x5[n, t]
                                         - 23 * m.x1[n, t] * m.x6[n, t] * m.u3[n])

m.ode7 = Constraint(m.n, m.tau,
                    rule=lambda m, n, t: m.dx7[n, t] == - m.q[n] * m.x7[n, t] + 46.0 * m.x1[n, t] * m.x6[n, t] * m.u3[n])

m.ode8 = Constraint(m.n, m.tau,
                    rule=lambda m, n, t: m.dx8[n, t] == 5.8 * (m.q[n] * m.x1[n, t] - m.u4[n]) - 3.7 * m.u1[n]
                    - 4.1 * m.u2[n] + m.q[n] * (23 * m.x4[n, t] + 11 * m.x5[n, t]
                    + 28 * m.x6[n, t] + 35 * m.x7[n, t]) - 5 * m.u3[n] ** 2 - 0.099)


def x1con(m, n):
    if n == m.n.first():
        return m.x1[n, 0] == 0.1883
    else:
        return m.x1[n, 0] == m.x1[n - 1, tauF]


m.x1connect = Constraint(m.n, rule=x1con)


def x2con(m, n):
    if n == m.n.first():
        return m.x2[n, 0] == 0.2507
    else:
        return m.x2[n, 0] == m.x2[n - 1, tauF]


m.x2connect = Constraint(m.n, rule=x2con)


def x3con(m, n):
    if n == m.n.first():
        return m.x3[n, 0] == 0.0467
    else:
        return m.x3[n, 0] == m.x3[n - 1, tauF]


m.x3connect = Constraint(m.n, rule=x3con)


def x4con(m, n):
    if n == m.n.first():
        return m.x4[n, 0] == 0.0899
    else:
        return m.x4[n, 0] == m.x4[n - 1, tauF]


m.x4connect = Constraint(m.n, rule=x4con)


def x5con(m, n):
    if n == m.n.first():
        return m.x5[n, 0] == 0.1804
    else:
        return m.x5[n, 0] == m.x5[n - 1, tauF]


m.x5connect = Constraint(m.n, rule=x5con)


def x6con(m, n):
    if n == m.n.first():
        return m.x6[n, 0] == 0.1394
    else:
        return m.x6[n, 0] == m.x6[n - 1, tauF]


m.x6connect = Constraint(m.n, rule=x6con)


def x7con(m, n):
    if n == m.n.first():
        return m.x7[n, 0] == 0.1046
    else:
        return m.x7[n, 0] == m.x7[n - 1, tauF]


m.x7connect = Constraint(m.n, rule=x7con)


def x8con(m, n):
    if n == m.n.first():
        return m.x8[n, 0] == 0.0000
    else:
        return m.x8[n, 0] == m.x8[n - 1, tauF]


m.x8connect = Constraint(m.n, rule=x8con)

m.obj = Objective(expr=m.x8[N, tauF], sense=maximize)
#discretizer = TransformationFactory('dae.finite_difference')
discretizer = TransformationFactory('dae.collocation')
discretizer.apply_to(m, wrt=m.tau, nfe=6)

solver = SolverFactory('ipopt')
solver.options['tol'] = 1e-6
solver.options['constr_viol_tol'] = 1e-6
results = solver.solve(m, tee=True)

u1 = [value(m.u1[n]) for n in m.n]

print(u1, value(m.x8[N, tauF]))
