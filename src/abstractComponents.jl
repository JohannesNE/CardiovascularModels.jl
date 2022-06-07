@connector function Con(;name)
    sts = @variables begin
        P(t) 
        Q(t), [connect = Flow]
    end
    ODESystem(Equation[], t, sts, []; name=name)
end

function CompliantCompartment(; name, PCargs...)

    @named pressurizedCompartment = PressurizedCompartment(;PCargs...)

    @named in = Con()
    @named out = Con()

    sts = @variables V(t) P(t) dV(t)

    # Flow is positive into the component.
    eqs = [
        D(V) ~ dV,
        dV ~ in.Q + out.Q,
        in.P ~ P,
        out.P ~ P
    ]
    
    pressurizedCompliantCompartment = extend(ODESystem(eqs, t, sts, []; name), 
                                            pressurizedCompartment)  
    compose(pressurizedCompliantCompartment, in, out)  
end

"""
This compartment can be thought of as a stiff tube with a fluid that has
mass. Inertia is the resistance to change in velocity of the fluid.

# Arguments
  - `Area`: [m^2] cross-sectional area of the compartment.
  - `V`: [m^3] Volume of the compartment. In the inertia compartment, Volume is a parameter (constant), not a state.
  - `ρ`: [kg/m^3] Density of the fluid.
  - `Q_start`: Initial condition for flow.
"""
function InertiaCompartment(; name, 
                               Area::Float64, 
                               V::Float64, 
                               ρ::Float64 = 1060., # kg/m^3.
                               Q_start::Float64 = 0.)

    @named volumelessComponent = VolumelessComponent()

    @unpack ΔP, Q = volumelessComponent

    @named in = Con()
    @named out = Con()
    
    ps = @parameters(
        V = V,
        Area = Area,
        ρ = ρ
    )    

    # Flow is positive into the component.
    eqs = [
        # Inductance equation.
        # L = ρ * V/Area^2 = rho * length / Area.
        D(Q) ~ ΔP / (ρ * V/Area^2)
    ]
    
    InertiaCompartment = extend(ODESystem(eqs, t, [], ps; name, defaults = [Q => Q_start]),
                                            volumelessComponent)  
    compose(InertiaCompartment, in, out)  
end

function VolumelessComponent(; name)
    @named in = Con()
    @named out = Con()

    sts = @variables Q(t) ΔP(t)

    # Flow is positive into the component.
    eqs = [
        0 ~ in.Q + out.Q,
        Q ~ in.Q,
        ΔP ~ in.P - out.P
    ]

    compose(ODESystem(eqs, t, sts, []; name), in, out)  
end

"""
# Arguments
  - `ext_pressure` and `tm_pressure`: If set to a number{Float64}, 
  this will fix the pressure to that value. 
  If set to `"free"`, no equation is made for controling
  the state (this should be done in an inheriting system (e.g. `CompliantCompartment()`)).
  
"""
function PressurizedCompartment(;name, ext_pressure = 0., tm_pressure = "free")

    sts = @variables P(t) P_ext(t) P_tm(t)

    eqs = [
        P ~ P_ext + P_tm
    ]

    if tm_pressure != "free"
        @assert typeof(tm_pressure) == Float64 

        push!(eqs, P_tm ~ tm_pressure)
    end

    if ext_pressure != "free"
        @assert typeof(ext_pressure) == Float64 

        push!(eqs, P_ext ~ ext_pressure)
    end

    ODESystem(eqs, t, sts, []; name)
end
PressurisedCompartment = PressurizedCompartment


