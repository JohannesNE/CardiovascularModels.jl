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
    Vd::Float64,
    ext_pressure=0.)

    @named compartment = VascularCompartment(;ext_pressure)
    @unpack V, P_tm = compartment

    ps = @parameters (Ees = Ees, Vd = Vd)

    eqs = [
        P_tm ~ max(0, Ees*(V-Vd))
    ]

    extend(ODESystem(eqs, t, [], ps; name), compartment)  
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
    # Arguments that are passed along to PressurisedCompartment
    ext_pressure = 0.,
    )

    @named compartment = VascularCompartment(;ext_pressure)
    @unpack V, P_tm = compartment

    ps = @parameters (
        Ees = Ees, 
        Vd = Vd, 
        V0 = V0, 
        λ = λ, 
        P0 = P0
    )

    sts = @variables (
        drv(t),
        Pes(t),
        Ped(t)
    )

    eqs = [
        Pes ~ Ees*(V-Vd),
        Ped ~ P0*(exp(λ*(V-V0))-1.0),
        P_tm ~ drv*Pes+Ped
    ]

    extend(ODESystem(eqs, t, [drv], ps; name), compartment)  
    
end

function Const_Pressure(; name, 
    # Parameters
    P::Float64)

    @named out = Con()

    ps = @parameters (P = P)

    eqs = [
        out.P ~ P
    ]

    compose(ODESystem(eqs, t, [], ps; name), out)  

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
