include("./const.jl")

module StructsMaxwell

using ..Constants
using ..Constants_Maxwell

export FDTD
export OpticsProperties
export CellType


Base.@kwdef mutable struct FDTD

    drude_region::Array{Bool,2} = zeros(Bool, Nx+1, Ny+1)

    Ez::Array{ComplexF64,2} = zeros(ComplexF64, Nx+1, Ny+1)
    Hx::Array{ComplexF64,2} = zeros(ComplexF64, Nx+1, Ny)
    Hy::Array{ComplexF64,2} = zeros(ComplexF64, Nx, Ny+1)
    
    σx::Array{ComplexF64,1} = zeros(ComplexF64, Nx)
    σy::Array{ComplexF64,1} = zeros(ComplexF64, Ny)
    
    ρₑy::Array{ComplexF64,1} = zeros(ComplexF64, PMLy*2)
    ρₘy::Array{ComplexF64,1} = zeros(ComplexF64, PMLy*2)
    ρₑx::Array{ComplexF64,1} = zeros(ComplexF64, PMLy*2)
    ρₘx::Array{ComplexF64,1} = zeros(ComplexF64, PMLy*2)
    
    σpₑy::Array{ComplexF64,1} = zeros(ComplexF64, PMLy*2)
    σpₘy::Array{ComplexF64,1} = zeros(ComplexF64, PMLy*2)
    κₑy::Array{ComplexF64,1} = zeros(ComplexF64, PMLy*2)
    κₘy::Array{ComplexF64,1} = zeros(ComplexF64, PMLy*2)
    αₑy::Vector{ComplexF64}  = zeros(ComplexF64, PMLy*2)
    αₘy::Vector{ComplexF64}  = zeros(ComplexF64, PMLy*2)
    
    σpₑx::Array{ComplexF64,1} = zeros(ComplexF64, PMLx*2)
    σpₘx::Array{ComplexF64,1} = zeros(ComplexF64, PMLx*2)
    κₑx::Array{ComplexF64,1} = zeros(ComplexF64, PMLx*2)
    κₘx::Array{ComplexF64,1} = zeros(ComplexF64, PMLx*2)
    αₑx::Vector{ComplexF64}  = zeros(ComplexF64, PMLx*2)
    αₘx::Vector{ComplexF64}  = zeros(ComplexF64, PMLx*2)
    
    κx::Vector{ComplexF64} = zeros(ComplexF64, Nx)
    κy::Vector{ComplexF64} = zeros(ComplexF64, Ny)
    αₑ::Vector{ComplexF64}  = zeros(ComplexF64, Nx)
    αy::Vector{ComplexF64} = zeros(ComplexF64, Ny)
    
    ΨEzy::Array{ComplexF64,2} = zeros(ComplexF64, Nx+1, PMLy*2)
    ΨEzx::Array{ComplexF64,2} = zeros(ComplexF64, PMLx*2, Ny+1)
    ΨHxy::Array{ComplexF64,2} = zeros(ComplexF64, Nx+1, PMLy*2)
    ΨHyx::Array{ComplexF64,2} = zeros(ComplexF64, PMLx*2, Ny+1)
    CΨEzy::Array{ComplexF64,2} = zeros(ComplexF64, Nx+1, PMLy*2)
    CΨEzx::Array{ComplexF64,2} = zeros(ComplexF64, PMLx*2, Ny+1)
    CΨHxy::ComplexF64 = 0 + 0im
    CΨHyx::ComplexF64 = 0 + 0im
    
    bEy::Array{ComplexF64,1} = zeros(ComplexF64, PMLy*2)
    aEy::Array{ComplexF64,1} = zeros(ComplexF64, PMLy*2)
    bMy::Array{ComplexF64,1} = zeros(ComplexF64, PMLy*2)
    aMy::Array{ComplexF64,1} = zeros(ComplexF64, PMLy*2)
    
    bEx::Array{ComplexF64,1} = zeros(ComplexF64, PMLx*2)
    aEx::Array{ComplexF64,1} = zeros(ComplexF64, PMLx*2)
    bMx::Array{ComplexF64,1} = zeros(ComplexF64, PMLx*2)
    aMx::Array{ComplexF64,1} = zeros(ComplexF64, PMLx*2)


    CHxEz::Array{ComplexF64,1} = zeros(Ny)  # Coefficient for Hx and Ez update
    CHyEz::Array{ComplexF64,1} = zeros(Nx)  # Coefficient for Hy and Ez update
    CEzHx::Array{ComplexF64,2} = zeros(Nx-1,Ny-1)  # Coefficient for Ez and Hx update
    CEzHy::Array{ComplexF64,2} = zeros(Nx-1,Ny-1)  # Coefficient for Ez and Hy update
    CEzJ::Array{ComplexF64,2} = zeros(Nx,Ny)  # Coefficient for Ez and J update

 

    ϵz::Array{ComplexF64,2} = ones(Nx, Ny)

    Jz::Array{ComplexF64,2} = zeros(ComplexF64, Nx+1, Ny+1)  # Current density (same size as Ez)
    aJ::Array{Float64,2} = zeros(Nx+1, Ny+1)                 # Coefficient for J update: (1 - γΔt/2)/(1 + γΔt/2)
    bJ::Array{Float64,2} = zeros(Nx+1, Ny+1)                 # Coefficient for J update: (ϵ₀ωₚ²Δt)/(1 + γΔt/2)
    ωp::Array{Float64,2} = zeros(Nx+1, Ny+1)                 # Plasma frequency [rad/s] (spatially varying)
    γ::Array{Float64,2} = zeros(Nx+1, Ny+1)                  # Collision frequency [1/s] (spatially varying)
    ϵ∞::Array{Float64,2} = ones(Nx+1, Ny+1)                  # High-frequency permittivity (default = 1)



end

Base.@kwdef mutable struct OpticsProperties
    n::Array{Float64,2} = ones(Nx,Ny)
    k::Array{Float64,2} = zeros(Nx,Ny)
    R::Array{Float64,2} = zeros(Nx,Ny)
    ω::Array{Float64,2} = ω₀*ones(Nx,Ny)
    ωₚ::Array{Float64,2} = zeros(Nx,Ny)
    λ::Array{Float64,2} = λ₀*ones(Nx,Ny)
    ν::Array{Float64,2} = zeros(Nx,Ny)
    γ::Array{Float64,2} = zeros(Nx,Ny)
end

end