#left Ventricle
@named l_ventricle = Ventricle(Ees = 200e6,
Vd = 0., V0 = 0., λ = 33e3, P0 = 10.,
R_out = 6e6)

@named card_driver = Driver()

# Aorta
@named aorta = Vessel(Ees = 200e6, Vd = 100e-6, R_out = 50e6)

eqs_driver = [
l_ventricle.drv ~ card_driver.contraction
]

#@named vein = Const_Pressure(P = 1e3, R_out = 6e6, valve_out = true)
@named vein = Vessel(Ees = 5e6, Vd = 500e-6, R_out = 10e6, valve_out = true)

eqs_con = [
connect(vein.out, l_ventricle.in),
connect(l_ventricle.out, aorta.in),
connect(aorta.out, vein.in)
]

eqs_comb = [eqs_driver; eqs_con]

# ODE parameters
volume_start = [l_ventricle.V => 1e-4, 
aorta.V => 1.5e-4,
vein.V => 2e-3]
time_span = (0.0, 10.0)

@named hemo_sys = ODESystem(eqs_comb)

@named connected = compose(hemo_sys, l_ventricle, aorta, vein, card_driver)

problem = ODEProblem(structural_simplify(connected), volume_start, time_span, [])

sol = solve(problem, Tsit5(), dtmax = 0.01, reltol = 1e-6)