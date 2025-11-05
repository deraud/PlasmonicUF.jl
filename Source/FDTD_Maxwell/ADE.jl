module ADE

using ..StructsMaxwell
using ..Constants
using ..Constants_Maxwell

export InitializeADEY!
export InitializeADEX!

function InitializeADEY!(f::FDTD)
    σmaxy = σy*(ρy+1)/(150π*dy)
    Threads.@threads for i in 1:PMLy

        f.ρₑy[i] = ((PMLy-i+1)-0.25)/PMLy
        f.ρₑy[end-i+1] = (PMLy+1-i-0.25)/PMLy
        f.ρₘy[i] = ((PMLy-i+1)+0.25)/PMLy
        f.ρₘy[end-i+1] = (PMLy+1-i+0.25)/PMLy

        f.σpₑy[i] = σmaxy * f.ρₑy[i]^ρy
        f.σpₑy[end-i+1] = σmaxy * f.ρₑy[end-i+1]^ρy
        f.σpₘy[i] = (μₒ/ϵₒ)*σmaxy*f.ρₘy[i]^ρy
        f.σpₘy[end-i+1] = (μₒ/ϵₒ)*σmaxy*f.ρₘy[end-i+1]^ρy

        f.κₑy[i] = 1 + (κy-1)*f.ρₑy[i]^ρy
        f.κₑy[end-i+1] = 1 + (κy-1)*f.ρₑy[end-i+1]^ρy
        f.κₘy[i] = 1 + (κy-1)*f.ρₘy[i]^ρy
        f.κₘy[end-i+1] = 1 + (κy-1)*f.ρₘy[end-i+1]^ρy

        f.αₑy[i] = αy_min + (αy_max - αy_min) * (1-f.ρₑy[i])
        f.αₑy[end-i+1] = αy_min + (αy_max - αy_min) * (1-f.ρₑy[end-i+1])
        f.αₘy[i] = (μₒ/ϵₒ)*(αy_min + (αy_max - αy_min) * (1-f.ρₘy[i]))
        f.αₘy[end-i+1] = (μₒ/ϵₒ)*(αy_min + (αy_max - αy_min) * (1-f.ρₘy[end-i+1]))

        f.bEy[i] = ϵₒ*f.κₑy[i]/(ϵₒ*f.κₑy[i] + dt_maxwell*(f.σpₑy[i]+f.αₑy[i]*f.κₑy[i]))
        f.aEy[i] = -dt_maxwell*f.σpₑy[i]/(dy*f.κₑy[i]*(ϵₒ*f.κₑy[i]+dt_maxwell*(f.σpₑy[i]+f.αₑy[i]*f.κₑy[i])))
        f.bMy[i] = μₒ*f.κₘy[i]/(μₒ*f.κₘy[i] + dt_maxwell*(f.σpₘy[i]+f.αₘy[i]*f.κₘy[i]))
        f.aMy[i] = -dt_maxwell*f.σpₘy[i]/(dy*f.κₘy[i]*(μₒ*f.κₘy[i]+dt_maxwell*(f.σpₘy[i]+f.αₘy[i]*f.κₘy[i])))
    
        f.bEy[end-i+1] = ϵₒ*f.κₑy[end-i+1]/(ϵₒ*f.κₑy[end-i+1] + dt_maxwell*(f.σpₑy[end-i+1]+f.αₑy[end-i+1]*f.κₑy[end-i+1]))
        f.aEy[end-i+1] = -dt_maxwell*f.σpₑy[end-i+1]/(dy*f.κₑy[end-i+1]*(ϵₒ*f.κₑy[end-i+1]+dt_maxwell*(f.σpₑy[end-i+1]+f.αₑy[end-i+1]*f.κₑy[end-i+1])))
        f.bMy[end-i+1] = μₒ*f.κₘy[end-i+1]/(μₒ*f.κₘy[end-i+1] + dt_maxwell*(f.σpₘy[end-i+1]+f.αₘy[end-i+1]*f.κₘy[end-i+1]))
        f.aMy[end-i+1] = -dt_maxwell*f.σpₘy[end-i+1]/(dy*f.κₘy[end-i+1]*(μₒ*f.κₘy[end-i+1]+dt_maxwell*(f.σpₘy[end-i+1]+f.αₘy[end-i+1]*f.κₘy[end-i+1])))
    end
    Threads.@threads for i in 1:Nx+1
        for j in 1:(PMLy*2)
            f.CΨEzy[i,j] = f.CEzHx[1,1]*dy
        end
    end
    f.CΨHxy = f.CHxEz[1,1]*dy
    Threads.@threads for i in 1:Nx-1
        for j in 1:PMLy
            f.CEzHx[i,j] = f.CEzHx[i,j]/f.κₑy[j]
            f.CEzHx[i,end-j+1] = f.CEzHx[i,end-j+1]/f.κₑy[end-j+1]
        end
    end

    Threads.@threads for j in 1:PMLy
        f.CHxEz[j] = f.CHxEz[j]/f.κₘy[j]
        f.CHxEz[Ny-j+1] = f.CHxEz[Ny-j+1]/f.κₘy[end-j+1]
    end
end

function InitializeADEX!(f::FDTD)
    σmaxx = σx * (ρx + 1) / (150π * dx)
    Threads.@threads for i in 1:PMLx
        f.ρₑx[i] = ((PMLx - i + 1) - 0.25) / PMLx
        f.ρₑx[end - i + 1] = ((PMLx + 1 - i) - 0.25) / PMLx
        f.ρₘx[i] = ((PMLx - i + 1) + 0.25) / PMLx
        f.ρₘx[end - i + 1] = ((PMLx + 1 - i) + 0.25) / PMLx

        f.σpₑx[i] = σmaxx * f.ρₑx[i]^ρx
        f.σpₑx[end - i + 1] = σmaxx * f.ρₑx[end - i + 1]^ρx
        f.σpₘx[i] = (μₒ / ϵₒ) * σmaxx * f.ρₘx[i]^ρx
        f.σpₘx[end - i + 1] = (μₒ / ϵₒ) * σmaxx * f.ρₘx[end - i + 1]^ρx

        f.κₑx[i] = 1 + (κx - 1) * f.ρₑx[i]^ρx
        f.κₑx[end - i + 1] = 1 + (κx - 1) * f.ρₑx[end - i + 1]^ρx
        f.κₘx[i] = 1 + (κx - 1) * f.ρₘx[i]^ρx
        f.κₘx[end - i + 1] = 1 + (κx - 1) * f.ρₘx[end - i + 1]^ρx

        f.αₑx[i] = αx_min + (αx_max - αx_min) * (1 - f.ρₑx[i])
        f.αₑx[end - i + 1] = αx_min + (αx_max - αx_min) * (1 - f.ρₑx[end - i + 1])
        f.αₘx[i] = (μₒ / ϵₒ) * (αx_min + (αx_max - αx_min) * (1 - f.ρₘx[i]))
        f.αₘx[end - i + 1] = (μₒ / ϵₒ) * (αx_min + (αx_max - αx_min) * (1 - f.ρₘx[end - i + 1]))

        f.bEx[i] = ϵₒ * f.κₑx[i] / (ϵₒ * f.κₑx[i] + dt_maxwell * (f.σpₑx[i] + f.αₑx[i] * f.κₑx[i]))
        f.aEx[i] = -dt_maxwell * f.σpₑx[i] / (dx * f.κₑx[i] * (ϵₒ * f.κₑx[i] + dt_maxwell * (f.σpₑx[i] + f.αₑx[i] * f.κₑx[i])))
        f.bMx[i] = μₒ * f.κₘx[i] / (μₒ * f.κₘx[i] + dt_maxwell * (f.σpₘx[i] + f.αₘx[i] * f.κₘx[i]))
        f.aMx[i] = -dt_maxwell * f.σpₘx[i] / (dx * f.κₘx[i] * (μₒ * f.κₘx[i] + dt_maxwell * (f.σpₘx[i] + f.αₘx[i] * f.κₘx[i])))

        f.bEx[end - i + 1] = ϵₒ * f.κₑx[end - i + 1] / (ϵₒ * f.κₑx[end - i + 1] + dt_maxwell * (f.σpₑx[end - i + 1] + f.αₑx[end - i + 1] * f.κₑx[end - i + 1]))
        f.aEx[end - i + 1] = -dt_maxwell * f.σpₑx[end - i + 1] / (dx * f.κₑx[end - i + 1] * (ϵₒ * f.κₑx[end - i + 1] + dt_maxwell * (f.σpₑx[end - i + 1] + f.αₑx[end - i + 1] * f.κₑx[end - i + 1])))
        f.bMx[end - i + 1] = μₒ * f.κₘx[end - i + 1] / (μₒ * f.κₘx[end - i + 1] + dt_maxwell * (f.σpₘx[end - i + 1] + f.αₘx[end - i + 1] * f.κₘx[end - i + 1]))
        f.aMx[end - i + 1] = -dt_maxwell * f.σpₘx[end - i + 1] / (dx * f.κₘx[end - i + 1] * (μₒ * f.κₘx[end - i + 1] + dt_maxwell * (f.σpₘx[end - i + 1] + f.αₘx[end - i + 1] * f.κₘx[end - i + 1])))
    end

    Threads.@threads for j in 1:Ny+1
        for i in 1:(PMLx*2)
            f.CΨEzx[i,j] = f.CEzHy[1,1] * dx
        end
    end
    f.CΨHyx = f.CHyEz[1,1] * dx

    Threads.@threads for j in 1:Ny-1
        for i in 1:PMLx
            f.CEzHy[i,j] = f.CEzHy[i,j] / f.κₑx[i]
            f.CEzHy[end - i + 1, j] = f.CEzHy[end - i + 1, j] / f.κₑx[end - i + 1]
        end
    end

    Threads.@threads for i in 1:PMLx
        f.CHyEz[i] = f.CHyEz[i] / f.κₘx[i]
        f.CHyEz[Nx - i + 1] = f.CHyEz[Nx - i + 1] / f.κₘx[end - i + 1]
    end
end

end