module Constants_Maxwell

using ..Constants

#include(joinpath(@__DIR__, "..", "const.jl"))
include("C:\\Users\\62812\\OneDrive\\Desktop\\Files\\ITB\\S2\\Nanobubble2\\const.jl")

export λ₀, c0, ω₀, CL, Nt, Nx, Ny, dx, dy, SL, DL, dt_calc, dt_maxwell,
       ϵₒ, μₒ, ϵᵣ,
       PMLx, PMLy, ρx, ρy, σx, σy, κx, κy,
       αx_min, αy_min, αx_max, αy_max,
       w0, sourcey, z0, z, zR, w_z, R_z, w, amplitude, psi_z, k0, tₚ

function floor_to_1_or_5(x::Real)
    if x == 0
        return 0
    end

    order = floor(Int, log10(abs(x)))      # order of magnitude
    base = 10.0^order                      # base scale
    scaled = x / base                      # scale x to [1,10)

    if scaled >= 5
        floorval = 5
    else
        floorval = 1
    end

    return sign(x) * floorval * base
end

const λ₀::Float64 = 800e-9
const c0::Float64 = 299792458.0  # Speed of light in vacuum (m/s)
const ω₀::Float64 = 2π*c0/λ₀
const CL::Float64 = 0.99  # Courant number for stability

#const Nt::Int = 1000  # Number of time steps
#const Nx::Int = 600  # Number of grid points in x-direction
#const Ny::Int = 1200  # Number of grid points in y-direction
const dx::Float64 = 40e-9#5.1944e-08  # Spatial step size in x-direction (m)   
const dy::Float64 = 1e-9#5.1944e-08  # Spatial step size in y-direction (m)

const SL::Float64 = 0.15  # Source Location 
const FL::Float64 = 0.65  # Focal Location
const DL::Float64 = 0.0  # Dielectric Location (0 = full dielectric, 1 = no dielectric)

#const dt_calc::Float64 = CL/(c0*sqrt(dx^-2 + dy^-2))
#const dt::Float64 = floor_to_1_or_5(dt_calc)
const dt::Float64 = 1e-17
const dt_calc::Float64 = 0.99/(c0*sqrt(dx^-2 + dy^-2))
const dtscale::Int = ceil(dt/dt_calc)
const dt_maxwell::Float64 = dt/dtscale # Time step size for FDTD (s)
#const dt_maxwell::Float64 = 1e-17

#const dt::Float64 = 1e-17
#const dt::Float64 = 0.99/(c0*sqrt(dx^-2 + dy^-2))

const ϵₒ::Float64 = 8.854*10^-12 # Permittivity of Vacuum [farad/meter]
const μₒ::Float64 = 4*pi*10^-7  # Permeability of Vacuum [henry/meter]
const k::Float64 = 2  # Extinction coefficient
const n::Float64 = 2  # Refractive index
const ϵᵣ::ComplexF64 = 1#(n + k*im)^2  # Relative Permittivity of the medium (e.g., steel) [farad/meter]

## PML Constants
const PMLx::Int = 20  # PML width in x-direction
const PMLy::Int = 20  # PML width in y-direction
const ρx::Float64 = 3.65  # Polynomial order in x-direction
const ρy::Float64 = 3.65  # Polynomial order in y-direction
const σx::Float64 = 1  # Conductivity in x-direction
const σy::Float64 = 1  # Conductivity in y-direction
const κx::Float64 = 2 # kappa in x-direction
const κy::Float64 = 2 # kappa in y-direction
const αx_min::Float64 = 0  # alpha min in x-direction
const αy_min::Float64 = 0  # alpha min in y-direction
const αx_max::Float64 = 4e-5  # alpha max in x-direction
const αy_max::Float64 = 4e-5  # alpha max in y-direction

## Laser Source
const tₚ::Float64 = 100e-15
const w0::Float64 = 1e-6  # Beam waist (m)
const k0::Float64 = ω₀ / c0  # Wave number
const sourcey::Int = PMLy+1 #trunc(Int,SL*Ny) # Source position 
const z0::Float64 = (120)*dy#FL*Ny*dy # Focal Point 
const z::Float64 = sourcey*dy - z0 # Source to Focal Point distance
const zR::Float64 = π*w0^2/λ₀  # Rayleigh range
const w_z::Float64 = w0 * sqrt(1 + (z/zR)^2)  # Beam radius at distance z
const R_z::Float64 = z * (1 + (zR / z)^2) # Radius of Curvature
const w::Float64 = w0 * sqrt(1 + (z / zR)^2) # Beam spot size
const amplitude::Float64 = w0/w_z
const psi_z::Float64 = atan(z / zR) # Gouy phase

end