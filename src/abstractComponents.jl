@connector function Con(;name)
    sts = @variables begin
        P(t) 
        Q(t), [connect = Flow]
    end
    ODESystem(Equation[], t, sts, []; name=name)
end

function VascularCompartment(; name)

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

function PressurizedCompartment(;name, P_ext = 0., P_m)

    sts = @variables P(t) P_ext(t) = P_ext P_m(t) = P_m

    eqs = [
        P ~ P_ext + P_m
    ]

    ODESystem(eqs, t, sts, []; name)
end
PressurisedCompartment = PressurizedCompartment


