module FDTDMaxwell

include(joinpath(@__DIR__, "..", "const.jl"))
using ..Constants
include("././const.jl")
include("././struct.jl")
include("././fdtd.jl")
include("././source.jl")
include("././ADE.jl")


using .Functions
using .ADE
using .Source
using .Constants_Maxwell
    
end