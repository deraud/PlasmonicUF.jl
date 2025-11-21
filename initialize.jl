module Initialization

include("./const.jl")
using .Constants

export InitializePhase!

function InitializePhase!(phase, R::Float64)

    center_x = Nx / 2.0
    center_y = Ny / 2.0

    for i in 1:Nx
        for j in 1:Ny

            dx = i - center_x
            dy = j - center_y
            distance = sqrt(dx^2 + dy^2)
            
            if distance <= R
                phase.cell_type[i,j] = CellType.SOLID
            else
                phase.cell_type[i,j] = CellType.GAS
            end
        end
    end 
end

end