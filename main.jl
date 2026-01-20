using Plots
using StaticArrays
using ProgressBars
using LinearAlgebra
using Statistics
# Optional for benchmarking (comment out if not installed)
# using BenchmarkTools

include("./const.jl")
include("./struct.jl")
include("./MRT3/MRT3.jl")
##include("./Source/FDTD_Maxwell/FDTD_Maxwell.jl")
include("./initialize.jl")
include("./animate.jl")


using .Constants
using .Struct
using .Initialization
using .MRT3
##using .FDTDMaxwell
using .Animate

"""
    run_sim!(lb, Nt, NtScale; show_progress=true, record=true, RhoEvolution=nothing)

Runs the MRT3 simulation in a function (faster than global scope).
- If `record==true`, returns a 3D array `RhoEvolution` of size (Nx, Ny, fld(Nt,NtScale)).
- Pass a preallocated `RhoEvolution` to reuse memory (must match expected size).
- Set `show_progress=false` when benchmarking to avoid I/O overhead.
"""
function run_sim!(lb::MRT3.Struct_MRT3.Lattice, Nt::Integer, NtScale::Integer;
                  show_progress::Bool = true,
                  record::Bool = true,
                  RhoEvolution::Union{Nothing,AbstractArray} = nothing)

    @assert Nt > 0 && NtScale > 0
    nframes = fld(Nt, NtScale)

    Nx_local = Constants.Nx
    Ny_local = Constants.Ny

    if record
        if RhoEvolution === nothing
            RhoEvolution = zeros(eltype(lb.ρ), Nx_local, Ny_local, nframes)
        else
            sx, sy, sz = size(RhoEvolution)
            @assert sx == Nx_local && sy == Ny_local "RhoEvolution first two dims must be ($Nx_local,$Ny_local), got ($sx,$sy)"
            @assert sz >= nframes "RhoEvolution has only $sz frames, needs at least $nframes"
        end
    end

    itr = show_progress ? ProgressBars.ProgressBar(1:Nt) : (1:Nt)

    for t in itr
        MRT3.EOS_SC!(lb)
        MRT3.calculate_force!(lb)
        MRT3.mrt_collision3!(lb)
        MRT3.streaming!(lb)
        MRT3.calculate_macro!(lb)

        if record && (t % NtScale == 0)
            frame = fld(t, NtScale)
            @inbounds @views RhoEvolution[:, :, frame] .= lb.ρ
        end
    end

    return record ? RhoEvolution : nothing
end


# ---------------------------
# Setup
# ---------------------------
lb = MRT3.Struct_MRT3.Lattice()
phase = Struct.Phase()

MRT3.initialize!(lb)

Initialization.InitializePhase!(phase, 20.0)


nframes = fld(Nt, NtScale)
RhoEvolution = zeros(eltype(lb.ρ), Nx, Ny, nframes)


"""
    Ultrashort Laser with Finite Difference Time Domain (FDTD)
"""

"""
    Fluid Dynamics with MRT LBM
"""
run_sim!(lb, NtScale, NtScale; show_progress=false, record=true, RhoEvolution=RhoEvolution)
RhoEvolution = run_sim!(lb, Nt, NtScale; show_progress=true, record=true, RhoEvolution=RhoEvolution)



println("Simulation Completed")

# Animate
#frames = 1:1:nframes
#animateDensityEvolution(RhoEvolution, frames)

# ---------------------------
# Optional: Benchmark
# ---------------------------
# @btime run_sim!($lb, $Nt, $NtScale; show_progress=false, record=true, RhoEvolution=$RhoEvolution);
