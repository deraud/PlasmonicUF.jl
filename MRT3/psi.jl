module Psi

using ..Constants
using ..Struct
using ..Struct_MRT3
using ..Constants_MRT3

export EOS_SC!, calculate_force!

function EOS_SC!(lb::Lattice)
    Threads.@threads for i in 1:Nx 
        for j in 1:Ny
            lb.Ψ[i,j] = ρ0 * (1.0 - exp(-lb.ρ[i,j]/ρ0))
        end
    end
end

function calculate_force!(lb::Lattice)
    Threads.@threads for i in 1:Nx
        for j in 1:Ny
            fx, fy = 0.0, 0.0
            for k in 1:Q
                ii = mod1(i + ex[k], Nx)
                jj = mod1(j + ey[k], Ny)
                fx += w[k] * ex[k] * lb.Ψ[ii, jj]
                fy += w[k] * ey[k] * lb.Ψ[ii, jj]
            end
            lb.Fx[i, j] = -G * lb.Ψ[i, j] * fx
            lb.Fy[i, j] = -G * lb.Ψ[i, j] * fy
        end
    end
end

end