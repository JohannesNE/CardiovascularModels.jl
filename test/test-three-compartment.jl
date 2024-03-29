using OrdinaryDiffEq, ModelingToolkit, CardiovascularModels

@variables t
D = Differential(t)

#left Ventricle
@named l_ventricle = Ventricle(Ees = 200e6,
    Vd = 0., V0 = 0., λ = 33e3, P0 = 10.)

@named card_driver = automaticDriver()

@named aortic_valve = Valve(R = 6e6)

# Aorta
@named aorta = ElasticVessel(Ees = 200e6, Vd = 100e-6)

@named systemic_resistance = Resistor(R = 50e6)

eqs_driver = [
    l_ventricle.drv ~ card_driver.contraction
]

@named vein = ElasticVessel(Ees = 5e6, Vd = 500e-6)
@named mitral_valve = Valve(R = 10e6)

eqs_con = serial_connect(vein, 
    mitral_valve, 
    l_ventricle,
    aortic_valve,
    aorta,
    systemic_resistance,
    vein)

eqs_comb = [eqs_driver; eqs_con]

# ODE parameters
volume_start = [l_ventricle.V => 1e-4, 
    aorta.V => 1.5e-4,
    vein.V => 2e-3]
time_span = (0.0, 10.0)

@named hemo_sys = ODESystem(eqs_comb)

@named connected = compose(hemo_sys, l_ventricle, card_driver, 
                            aortic_valve, aorta, systemic_resistance, vein, mitral_valve)

problem = ODEProblem(structural_simplify(connected), volume_start, time_span, [])

sol = solve(problem, Tsit5(), dtmax = 0.01, reltol = 1e-6)

# Plots.plot(sol)