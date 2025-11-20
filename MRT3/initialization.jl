module Initialization_MRT3

using ..Constants
using ..Struct
using ..Struct_MRT3
using ..Constants_MRT3

export initialize!

function initialize!(lb::Lattice)
    cx, cy = Nx÷2, Ny÷2
    radius = min(Nx, Ny) ÷ 6
    
    ρ_liquid = 7.2 * ρ0  # High density liquid phase
    ρ_vapor = 0.2 * ρ0   # Low density vapor phase
    interface_width = 3.0  # Width of interface transition
    
    for i in 1:Nx, j in 1:Ny
        dist = sqrt((i-cx)^2 + (j-cy)^2)
        
        # Smooth interface using hyperbolic tangent
        lb.ρ[i, j] = ρ_vapor + 0.5 * (ρ_liquid - ρ_vapor) * 
                  (1.0 - tanh((dist - radius) / interface_width))
        
        # Initialize distribution functions to equilibrium (zero velocity)
        lb.u[i, j] = 0.0
        lb.v[i, j] = 0.0
        for k in 1:Q
            lb.f[i, j, k] = w[k] * lb.ρ[i, j]
        end
    end
end

end