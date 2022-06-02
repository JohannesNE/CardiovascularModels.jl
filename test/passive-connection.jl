using ModelingToolkit, OrdinaryDiffEq, CardiovascularModels

@variables t
D = Differential(t)

@named vessel1 = ElasticVessel(Ees = 5e6, Vd = 0.0, ext_pressure = 1e2)
@named resistor = Resistor(R = 10e6)
@named vessel2 = ElasticVessel(Ees = 1e7, Vd = 0.0)

eqs_con = [
connect(vessel1.out, resistor.in),
connect(resistor.out, vessel2.in),
vessel1.in.Q ~ 0,
vessel2.out.Q ~ 0
]

volume_start = [
    vessel1.V => 1e-4, 
    vessel2.V => 2e-4
]

time_span = (0.0, 5.0)

@named connect_sys = ODESystem(eqs_con, t)
@named test1 = compose(connect_sys, vessel1, vessel2, resistor)
problem = ODEProblem(structural_simplify(test1), volume_start, time_span, [])

sol = solve(problem, Tsit5(), dtmax = 0.01, reltol = 1e-6)

# Plots.plot(sol)

# Test Inertia

@named inertia = InertiaCompartment(Area = 0.4e-3, V = 0.2e-3)

eqs_con_inertia = [
connect(vessel1.out, inertia.in),
connect(inertia.out, vessel2.in),
vessel1.in.Q ~ 0,
vessel2.out.Q ~ 0
]

@named connect_sys_inertia = ODESystem(eqs_con_inertia, t)
@named test_inertia = compose(connect_sys_inertia, vessel1, vessel2, inertia)
problem_inertia = ODEProblem(structural_simplify(test_inertia), volume_start, time_span, [])

sol_inertia = solve(problem_inertia, Tsit5(), dtmax = 0.01, reltol = 1e-6)

# Plots.plot(sol_inertia)
