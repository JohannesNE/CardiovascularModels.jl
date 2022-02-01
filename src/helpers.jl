
function serial_connect(syss::Vector)
    n_sys = length(syss)
    @assert n_sys > 1 "Connect at least two components"
    
    connections = [
        connect(syss[i].out, syss[i+1].in) 
        for i in 1:(n_sys-1)
    ]

    return connections
end

serial_connect(syss...) = serial_connect([syss...])