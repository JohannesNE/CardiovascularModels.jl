# Test Inertia

using ModelingToolkit, OrdinaryDiffEq, CardiovascularModels

@variables t
D = Differential(t)

@named vessel1 = DampedElasticVessel(Ees = 5e6, Vd = 0.0, R_damp = 1e6)
@named vessel2 = DampedElasticVessel(Ees = 1e7, Vd = 0.0, R_damp = 1e6)

@named inertia = InertiaCompartment(Area = 0.4e-3, V = 0.2e-3)


eqs_con_inertia = [
    connect(vessel1.out, inertia.in),
    connect(inertia.out, vessel2.in),
    vessel1.in.Q ~ 0,
    vessel2.out.Q ~ 0
    ]
    
time_span = (0.0, 5.0)
@named connect_sys_inertia = ODESystem(eqs_con_inertia, t)
@named test_inertia = compose(connect_sys_inertia, vessel1, vessel2, inertia)
problem_inertia = ODEProblem(structural_simplify(test_inertia), volume_start, time_span, [])

sol_inertia = solve(problem_inertia, Rodas5(), dtmax = 0.01, reltol = 1e-6)

# Plots.plot(sol_inertia, vars = [vessel1.V, vessel2.V])
