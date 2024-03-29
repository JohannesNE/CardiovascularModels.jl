# Simple 3-compartment model

A simple model with a ventricle, an artery (aorta) and a vein, connected with some resistance between them.

```julia
using OrdinaryDiffEq, ModelingToolkit, CardiovascularModels, Plots

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

# Connect systems in series.
# vein.out -> mitral_valve.in
# mitral_valve.out -> l_ventricle.in
# etc.
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
```


```julia
pascal2mmhg(t, pascal) = t, pascal * 0.00750062
m32ml(t, m3) = t, m3 * 1e6
```

```
m32ml (generic function with 1 method)
```





## Volume plots [ml]
```julia
plot_v1 = plot(sol, vars=[(m32ml, 0, l_ventricle.V), (m32ml, 0, aorta.V)], tspan =  (1,6));
plot_v2 = plot(sol, vars=[(m32ml, 0, vein.V)], tspan =  (1,6));
plot(plot_v1, plot_v2, layout = (2,1), ylabel = "ml")
```

![](figures/three_compartment_model_3_1.png)



## Pressure plots [mmHg]
```julia
plot_p1 = plot(sol, vars=[(pascal2mmhg, 0,l_ventricle.P), (pascal2mmhg, 0,aorta.P)], tspan =  (1,6));
plot_p2 = plot(sol, vars=[(pascal2mmhg, 0, vein.P)], tspan =  (1,6));
plot(plot_p1, plot_p2, layout = (2,1), ylabel = "mmHg")
```

![](figures/three_compartment_model_4_1.png)
