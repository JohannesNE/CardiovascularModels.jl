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

    @variables V(t) P(t)

    # Flow is positive into the component.
    eqs = [
        D(V) ~ in.Q + out.Q,
        in.P ~ P,
        in.P ~ out.P
    ]
    
    pressurizedCompliantCompartment = extend(ODESystem(eqs, t, [V, P], []; name), 
                                            pressurizedCompartment)  
    compose(pressurizedCompliantCompartment, in, out)  
end

function InertialCompartment(; name, PCargs...)

    @named pressurizedCompartment = PressurizedCompartment(;PCargs...)

    @named in = Con()
    @named out = Con()

    @variables V(t) P(t) 

    # Flow is positive into the component.
    eqs = [
        D(V) ~ in.Q + out.Q,
        ΔP ~ in.P - out.P
    ]
    
    pressurizedCompliantCompartment = extend(ODESystem(eqs, t, [V, P], []; name), 
                                            pressurizedCompartment)  
    compose(pressurizedCompliantCompartment, in, out)  
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


