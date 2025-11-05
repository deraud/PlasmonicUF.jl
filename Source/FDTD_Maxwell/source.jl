module Source

using ..StructsMaxwell
using ..Constants
using ..Constants_Maxwell

export CalculateGaussianSource!

function CalculateGaussianSource!(f::FDTD, t::Int)
    for i in PMLx+1:Nx-PMLx-1
        r2 = (i*dx-trunc(Int, Nx/2)*dx)^2
        if z != 0
            phase = k0 * (z + r2 / (2 * R_z)) - psi_z
        else
            phase = -psi_z
        end
        spatial = exp(-r2/w_z^2)
        temporal = sin(ω₀*t*dt_maxwell - phase)
        amplitude = w0/w_z
        mode_lock = exp(-(t*dt_maxwell - 2tₚ)^2 / tₚ^2)
        phys = exp( -2π*z0/λ₀*im - π*(i*dx-trunc(Int, Nx/2)*dx)^2*im/(2*R_z*λ₀) + psi_z*im)
        f.Ez[i,sourcey] += amplitude * spatial * temporal * phys * mode_lock
    end
end

end