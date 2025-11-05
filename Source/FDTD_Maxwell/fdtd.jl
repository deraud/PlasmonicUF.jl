module Functions

using ..StructsMaxwell
using ..Constants
using ..Constants_Maxwell

export InitializeUpdate!
export UpdateMagnetic!
export UpdateElectric!

function InitializeUpdate!(f::FDTD, phase)
    for i in 1:Nx
        for j in 1:Ny
            if (phase.cell_type[i, j] & CellType.SOLID) == 0
                continue
            end
            f.drude_region[i,Ny+1-j] = true
            f.ωp[i,Ny+1-j] = 1.37e16
            f.γ[i,Ny+1-j] = 4.04e13
            
        end
    end
    #i0 = div(Nx + 1, 2)   # x-center of crater
    #j0 = 300              # surface of the metal
    #R = 200                # radius of the hemisphere in grid units
    #for i in 1:Nx+1
    #    for j in 300:Ny+1
    #        if (i - i0)^2 + (j - j0)^2 <= R^2
    #            f.drude_region[i,j] = false
    #            f.ωp[i,j] = 0.0
    #            f.γ[i,j] = 0.0
    #        end
    #    end
    #end
    # --- Update coefficients for H fields (unchanged) ---
    Threads.@threads for i in 1:Nx
        f.CHyEz[i] = 2*dt_maxwell/(2*μₒ*dx)
    end
    Threads.@threads for j in 1:Ny
        f.CHxEz[j] = -2*dt_maxwell/(2*μₒ*dy)
    end

    # --- Initialize material regions ---
    Threads.@threads for i in PMLx:Nx-PMLx
        if DL != 0
            # Set background permittivity (non-Drude region)
            for j in trunc(Int,Ny*DL):Ny
                f.ϵz[i,j] = ϵᵣ
                # Default Drude parameters to 0 (non-Drude)
                f.ωp[i,j] = 0.0
                f.γ[i,j] = 0.0
                f.ϵ∞[i,j] = 1.0
            end
        end
    end

    # --- Combined FDTD + Drude coefficients ---
    Threads.@threads for i in 1:Nx-1
        for j in 1:Ny-1
            # For Drude regions: Use ϵ∞ instead of ϵz
            effective_ϵ = (f.ωp[i,j] > 0) ? f.ϵ∞[i,j] : f.ϵz[i,j]
            
            f.CEzHx[i,j] = -2*dt_maxwell/(2*ϵₒ*effective_ϵ*dy)
            f.CEzHy[i,j] = 2*dt_maxwell/(2*ϵₒ*effective_ϵ*dx)
            
            # Initialize Drude coefficients if in Drude region
            if f.ωp[i,j] > 0
                f.aJ[i,j] = (1 - f.γ[i,j] * dt_maxwell / 2) / (1 + f.γ[i,j] * dt_maxwell / 2)
                f.bJ[i,j] = (ϵₒ * f.ωp[i,j]^2 * dt_maxwell) / (1 + f.γ[i,j] * dt_maxwell / 2)
                f.CEzJ[i,j] = -dt_maxwell / (ϵₒ * f.ϵ∞[i,j])  # Negative sign for update
            end
        end
    end
end

function UpdateMagnetic!(f::FDTD)
    Threads.@threads for i in 1:Nx+1
        for j in 1:Ny
            f.Hx[i,j] += f.CHxEz[j] * (f.Ez[i,j+1] - f.Ez[i,j])
        end
    end
    
    Threads.@threads for i in 1:Nx
        for j in 1:Ny+1
            f.Hy[i,j] += f.CHyEz[i] * (f.Ez[i+1,j] - f.Ez[i,j])
        end
    end

    Threads.@threads for i in 1:(PMLx)
        for j in 1:Ny+1
            f.ΨHyx[i,j] = f.bMx[i] * f.ΨHyx[i,j] + f.aMx[i] * (f.Ez[i+1,j] - f.Ez[i,j])
            f.ΨHyx[end-i,j] = f.bMx[end-i] * f.ΨHyx[end-i] + f.aMx[end-i] * (f.Ez[end-i+1,j] - f.Ez[end-i,j])
        end
    end
    
    Threads.@threads for i in 1:(PMLx)
        for j in 1:Ny
            f.Hy[i,j] += f.CΨHyx * f.ΨHyx[i,j]
            f.Hy[end-i+1,j] += f.CΨHyx * f.ΨHyx[end-i+1,j]
        end
    end
    
    Threads.@threads for i in 1:Nx
        for j in 1:PMLy
            f.ΨHxy[i,j] = f.bMy[j] * f.ΨHxy[i,j] + f.aMy[j] * (f.Ez[i,j+1] - f.Ez[i,j])
            f.ΨHxy[i,end-j+1] = f.bMy[end-j+1] * f.ΨHxy[i,end-j+1] + f.aMy[end-j+1] * (f.Ez[i,end-j+1] - f.Ez[i,end-j])
        end
    end
    
    Threads.@threads for i in 1:Nx
        for j in 1:(PMLy)
            f.Hx[i,j] += f.CΨHxy * f.ΨHxy[i,j]
            f.Hx[i,end-j+1] += f.CΨHxy * f.ΨHxy[i,end-j+1]
        end
    end
end

function UpdateElectric!(f::FDTD)

   # Update Jz first (semi-implicit Drude model)
    Threads.@threads for i in 1:size(f.Jz, 1)
        for j in 1:size(f.Jz, 2)
            if f.drude_region[i, j]
                f.Jz[i, j] = f.aJ[i, j] * f.Jz[i, j] + f.bJ[i, j] * f.Ez[i, j]
            end
        end
    end

    # Standard FDTD update for Ez (with Drude correction)
    Threads.@threads for i in 2:Nx
        for j in 2:Ny
            f.Ez[i, j] += f.CEzHy[i-1, j-1] * (f.Hy[i, j] - f.Hy[i-1, j]) +
                          f.CEzHx[i-1, j-1] * (f.Hx[i, j] - f.Hx[i, j-1])
            # Add Drude current contribution
            if f.drude_region[i, j]
                f.Ez[i, j] += f.CEzJ[i, j] * f.Jz[i, j]
            end
        end
    end



    Threads.@threads for i in 1:PMLx
        for j in 1:Ny+1
            f.ΨEzx[i, j] = f.bEx[i] * f.ΨEzx[i, j] + 
                            f.aEx[i] * (f.Hy[i+1, j] - f.Hy[i, j])

            ii = Nx - PMLx + i -1 # actual index in Hy and Ez
            f.ΨEzx[PMLx + i, j] = f.bEx[PMLx + i] * f.ΨEzx[PMLx + i, j] +
                                  f.aEx[PMLx + i] * (f.Hy[ii+1, j] - f.Hy[ii, j])
        end
    end
  
    Threads.@threads for i in 1:PMLx
        ii = Nx - PMLx + i
        for j in 1:Ny+1
            f.Ez[i+1, j] += f.CΨEzx[i, j] * f.ΨEzx[i, j]
            f.Ez[ii, j] += f.CΨEzx[PMLx + i, j] * f.ΨEzx[PMLx + i, j]
        end
    end

    Threads.@threads for i in 1:Nx+1
        for j in 1:PMLy
            f.ΨEzy[i, j] = f.bEy[j] * f.ΨEzy[i, j] + f.aEy[j] * (f.Hx[i, j+1] - f.Hx[i, j])
            jj = Ny - PMLy + j -1
            f.ΨEzy[i, PMLy + j] = f.bEy[PMLy + j] * f.ΨEzy[i, PMLy + j] +
                                  f.aEy[PMLy + j] * (f.Hx[i, jj+1] - f.Hx[i, jj])
        end
    end

    Threads.@threads for i in 1:Nx+1
        for j in 1:PMLy
            f.Ez[i, j+1] += f.CΨEzy[i, j] * f.ΨEzy[i, j]
            jj = Ny - PMLy + j
            f.Ez[i, jj] += f.CΨEzy[i, PMLy + j] * f.ΨEzy[i, PMLy + j]
        end
    end
end

end
