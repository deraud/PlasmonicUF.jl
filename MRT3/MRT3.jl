module MRT3

include(joinpath(@__DIR__, "..", "const.jl"))
include(joinpath(@__DIR__, "..", "struct.jl"))

include("././const.jl")
using .Constants_MRT3
include("././struct.jl")
using .Struct_MRT3

include("./initialization.jl")
using .Initialization_MRT3
include("./lb.jl")
using .Lb_MRT
include("./Psi.jl")
using .Psi

end