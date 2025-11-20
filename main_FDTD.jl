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
include("./animate.jl")
include("./Source/FDTD_Maxwell/FDTD_Maxwell.jl")

using .Constants
using .Struct
using .MRT3
using .FDTDMaxwell
using .Animate