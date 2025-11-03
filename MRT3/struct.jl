include("./const.jl")
include(joinpath(@__DIR__, "..", "const.jl"))

module Struct_MRT3

using ..Constants_MRT3
using ..Constants

export Lattice

Base.@kwdef mutable struct Lattice
    ρ::Array{Float64,2} = zeros(Float64, Nx, Ny)
    f::Array{Float64,3} = zeros(Float64, Nx, Ny, Q)
    Fm::Array{Float64,3} = zeros(Float64, Nx, Ny, Q)  # force term in moment space
    u::Array{Float64,2} = zeros(Float64, Nx, Ny)
    v::Array{Float64,2} = zeros(Float64, Nx, Ny)
    pEOS::Array{Float64,2} = zeros(Float64, Nx, Ny)   # pressure from equation of state
    T::Array{Float64,2} = fill(300, Nx, Ny)             # temperature field
    Ψ::Array{Float64,2} = zeros(Float64, Nx, Ny)    # interaction potential
    Fx::Array{Float64,2} = zeros(Float64, Nx, Ny)   # force in x direction
    Fy::Array{Float64,2} = zeros(Float64, Nx, Ny)   # force in y direction
end

end