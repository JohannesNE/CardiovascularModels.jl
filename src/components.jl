@connector function Con(;name)
    sts = @variables begin
        P(t) 
        Q(t), [connect = Flow]
    end
    ODESystem(Equation[], t, sts, []; name=name)
end

function Compartment(; name)

    @named in = Con()
    @named out = Con()

    sts = @variables V(t) P(t)

    # Flow is positive into the component.
    eqs = [
        D(V) ~ in.Q + out.Q,
        in.P ~ P,
        in.P ~ out.P
    ]

    compose(ODESystem(eqs, t, [V, P], []; name), in, out)  
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

function Resistor(; name, R::Float64)
    @named component = VolumelessComponent()
    @unpack Q, ΔP = component

    ps = @parameters (R = R)

    eqs = [
        Q ~ ΔP / R
    ]

    extend(ODESystem(eqs, t, [Q, ΔP], ps; name), component)
end

function Valve(; name, R::Float64)
    @named component = VolumelessComponent()
    @unpack Q, ΔP = component

    ps = @parameters (R = R)

    eqs = [
        Q ~ (ΔP / R) * (ΔP > 0.0)
    ]

    extend(ODESystem(eqs, t, [Q, ΔP], ps; name), component)
end

"""
Vessel model (e.g. artery or vein)

# Arguments
   - `Ees`: Elasticity.
   - `Vd`: Unstressed volume (volume with no pressure).
"""
function Vessel(; name, 
    # Parameters
    Ees::Float64, 
    Vd::Float64)

    @named compartment = Compartment()
    @unpack V, P, in, out = compartment

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
    P0::Float64 = 0.
    )

    @named compartment = Compartment()
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
        Ped ~ P0*(exp(λ*(V-V0))-1.0),
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
_e_card_rep(t, A, B, C, cycle_len = 1.0) = _e_card(t % cycle_len, A, B, C)
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
