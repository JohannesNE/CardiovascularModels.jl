@connector function Con(;name)
    sts = @variables begin
        P(t) 
        Q(t), [connect = Flow]
    end
    ODESystem(Equation[], t, sts, []; name=name)
end

function Compartment(; name, 
    # Parameters
    R_out::Float64,
    valve_out::Bool = false)

    @named in = Con()
    @named out = Con()

    ps = @parameters (R_out = R_out)

    sts = @variables V(t) P(t)

    # Flow is positive into the compartment.
    if valve_out
        eq_Q_out = out.Q ~ -((P - out.P) / R_out) * (P > out.P)
    else
        eq_Q_out = out.Q ~ -((P - out.P) / R_out)
    end

    eqs = [
        eq_Q_out,
        D(V) ~ in.Q + out.Q,
        in.P ~ P
    ]

    compose(ODESystem(eqs, t, [V, P], ps; name), in, out)  
end

"""
Vessel model (e.g. artery or vein)

# Arguments
   - `Ees`: Elasticity.
   - `Vd`: Unstressed volume (volume with no pressure).
   - `R_out`: Resistance at outflow orifice.
   - `valve_out::Bool`: Valve at outflow orifice.

"""
function Vessel(; name, 
    # Parameters
    Ees::Float64, 
    Vd::Float64, 
    R_out::Float64,
    valve_out::Bool = false)

    @named compartment = Compartment(R_out = R_out, valve_out = valve_out)
    @unpack V, P = compartment

    ps = @parameters (Ees = Ees, Vd = Vd)

    eqs = [
        P ~ max(0, Ees*(V-Vd))
    ]

    extend(ODESystem(eqs, t, [V, P], ps; name), compartment)  

end

"""
Ventricle model

# Arguments
   - `Ees`: Maximum elasticity (systole).
   - `Vd`: (Unstressed volume) Volume where not contraction can occour.
   - `V0`: Volume with no passive elasticitiy.
   - `λ`: ?. Something diastolic
   - `P0`: ?. Something diastolic Unstressed
"""
function Ventricle(; name, 
    # Parameters
    Ees::Float64, 
    Vd::Float64, 
    V0::Float64, 
    λ::Float64,
    P0::Float64 = 0.,
    R_out::Float64,
    valve_out::Bool = true
    )

    @named compartment = Compartment(R_out = R_out, valve_out = valve_out)
    @unpack V, P = compartment

    ps = @parameters (
        Ees = Ees, 
        Vd = Vd, 
        V0 = V0, 
        λ = λ, 
        P0 = P0
    )

    sts = @variables (
        V(t), 
        P(t),
        drv(t),
        Pes(t),
        Ped(t)
    )

    eqs = [
        Pes ~ Ees*(V-Vd),
        Ped ~ P0*(exp(λ*(V-V0))-1),
        P ~ drv*Pes+Ped
    ]
    


    extend(ODESystem(eqs, t, [V, P, drv], ps; name), compartment)  
    
end

function Const_Pressure(; name, 
    # Parameters
    P::Float64, 
    R_out::Float64,
    valve_out::Bool = false)

    @named in = Con()
    @named out = Con()

    ps = @parameters (P = P)

    if valve_out
        eq_Q_out = out.Q ~ ((P - out.P) / R_out) * (P > out.P)
    else
        eq_Q_out = out.Q ~ ((P - out.P) / R_out)
    end

    eqs = [
        eq_Q_out,
        in.P ~ P
    ]

    compose(ODESystem(eqs, t, [], ps; name), in, out)  

end

# Register cardiac driver
_e_card(t, A, B, C) = A*exp(-B*(t-C)^2)
_e_card_rep(t, A, B, C, cycle_len = 1) = _e_card(t % cycle_len, A, B, C)
@register _e_card_rep(t, A, B, C) 

function Driver(; name,
    A::Float64 = 1., 
    B::Float64 = 80., 
    C::Float64 = 0.27)

    ps = @parameters (A=A, B=B, C=C)

    sts = @variables contraction(t)

    eqs = [
        contraction ~ _e_card_rep(t, A, B, C)
    ]

    ODESystem(eqs, t, [contraction], ps; name)  
end
