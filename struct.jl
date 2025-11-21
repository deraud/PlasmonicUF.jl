module Struct

include("./const.jl")
using .Constants

export Phase

Base.@kwdef mutable struct Phase
    cell_type::Matrix{Int} = fill(CellType.GAS, Nx, Ny)
end

end