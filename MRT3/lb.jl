module Lb_MRT

using ..Constants
using ..Struct
using ..Struct_MRT3
using ..Constants_MRT3

export mrt_collision!, streaming!, calculate_macro!
function mrt_collision!(lb::Lattice)
    Threads.@threads for i in 1:Nx 
        for j in 1:Ny
            # Add force effect to velocity (half time step)
            ux_eff = lb.u[i, j] + 0.5 * lb.Fx[i, j] / lb.ρ[i, j]
            uy_eff = lb.v[i, j] + 0.5 * lb.Fy[i, j] / lb.ρ[i, j]
            
            # Transform to moment space
            m = zeros(Q)
            for k in 1:Q
                for l in 1:Q
                    m[k] += M[k, l] * lb.f[i, j, l]
                end
            end

            # Calculate equilibrium moments
            usqr = ux_eff^2 + uy_eff^2
            meq = zeros(Q)
            meq[1] = lb.ρ[i, j]
            meq[2] = -2*lb.ρ[i, j] + 3*lb.ρ[i, j]*usqr
            meq[3] = lb.ρ[i, j] - 3*lb.ρ[i, j]*usqr
            meq[4] = lb.ρ[i, j]*ux_eff
            meq[5] = -lb.ρ[i, j]*ux_eff
            meq[6] = lb.ρ[i, j]*uy_eff
            meq[7] = -lb.ρ[i, j]*uy_eff
            meq[8] = lb.ρ[i, j]*(ux_eff^2 - uy_eff^2)
            meq[9] = lb.ρ[i, j]*ux_eff*uy_eff

            # Forcing term in moment space (Guo et al. 2002)
            F_dot_u = lb.Fx[i, j]*ux_eff + lb.Fy[i, j]*uy_eff
            Fm = zeros(Q)
            Fm[1] = 0.0
            Fm[2] = 6.0 * F_dot_u
            Fm[3] = -6.0 * F_dot_u
            Fm[4] = lb.Fx[i, j]
            Fm[5] = -lb.Fx[i, j]
            Fm[6] = lb.Fy[i, j]
            Fm[7] = -lb.Fy[i, j]
            Fm[8] = 2.0 * (lb.Fx[i, j]*ux_eff - lb.Fy[i, j]*uy_eff)
            Fm[9] = lb.Fx[i, j]*uy_eff + lb.Fy[i, j]*ux_eff

            # Collision in moment space with forcing
            Fm_scaled = (1.0 .- 0.5 .* s) .* Fm
            m_post = m - s .* (m - meq) + Fm_scaled

            # Transform back to velocity space
            for k in 1:Q
                lb.f[i, j, k] = 0.0
                for l in 1:Q
                    lb.f[i, j, k] += Minv[k, l] * m_post[l]
                end
            end
        end
    end
end
# MRT collision operator
function mrt_collisiona!(lb::Lattice)
    for i in 1:Nx, j in 1:Ny
        # Add force effect to velocity
        ux_eff = lb.u[i, j] + 0.5 * lb.Fx[i, j] / lb.ρ[i, j]
        uy_eff = lb.v[i, j] + 0.5 * lb.Fy[i, j] / lb.ρ[i, j]
        
        # Transform to moment space
        m = zeros(Q)
        for k in 1:Q
            for l in 1:Q
                m[k] += M[k, l] * lb.f[i, j, l]
            end
        end
        
        # Calculate equilibrium moments
        usqr = ux_eff^2 + uy_eff^2
        meq = zeros(Q)
        meq[1] = lb.ρ[i, j]
        meq[2] = -2*lb.ρ[i, j] + 3*lb.ρ[i, j]*usqr
        meq[3] = lb.ρ[i, j] - 3*lb.ρ[i, j]*usqr
        meq[4] = lb.ρ[i, j]*ux_eff
        meq[5] = -lb.ρ[i, j]*ux_eff
        meq[6] = lb.ρ[i, j]*uy_eff
        meq[7] = -lb.ρ[i, j]*uy_eff
        meq[8] = lb.ρ[i, j]*(ux_eff^2 - uy_eff^2)
        meq[9] = lb.ρ[i, j]*ux_eff*uy_eff
        
        # Collision in moment space
        m_post = m - s .* (m - meq)
        
        # Transform back to velocity space
        for k in 1:Q
            lb.f[i, j, k] = 0.0
            for l in 1:Q
                lb.f[i, j, k] += Minv[k, l] * m_post[l]
            end
            
            # Add force term
            ci_dot_F = ex[k]*lb.Fx[i, j] + ey[k]*lb.Fy[i, j]
            ci_dot_u = ex[k]*ux_eff + ey[k]*uy_eff
            lb.f[i, j, k] += w[k] * (1 - 0.5*s[1]) * 
                          (3*ci_dot_F + 9*ci_dot_F*ci_dot_u - 3*(lb.Fx[i,j]*ux_eff + lb.Fy[i,j]*uy_eff))
        end
    end
end

function streaming!(lb::Lattice)
    f_new = similar(lb.f)
    Threads.@threads for i in 1:Nx 
        for j in 1:Ny 
            for k in 1:Q
                ii = mod1(i + ex[k], Nx)
                jj = mod1(j + ey[k], Ny)
                f_new[ii, jj, k] = lb.f[i, j, k]
            end
        end
    end
    lb.f .= f_new
end

function calculate_macro!(lb::Lattice)
    Threads.@threads for i in 1:Nx 
        for j in 1:Ny
            lb.ρ[i, j] = sum(lb.f[i, j, :])
            lb.u[i, j] = sum(ex .* lb.f[i, j, :]) / lb.ρ[i, j]
            lb.v[i, j] = sum(ey .* lb.f[i, j, :]) / lb.ρ[i, j]
        end
    end
end

end