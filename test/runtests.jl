using CardiovascularModels
using Test

using OrdinaryDiffEq, ModelingToolkit

@variables t
D = Differential(t)

@testset "3-compartment model" begin
    include("test-three-compartment.jl")
end

@testset "Passive connection of two vessels" begin
    include("passive-connection.jl")
end