module Animate

using Plots

export animateDensityEvolution

function animateDensityEvolution(RhoEvolution, frames::Union{Vector{Int}, StepRange{Int, Int}})
    max = maximum(RhoEvolution)
    min = minimum(RhoEvolution)
    anim = @animate for t in frames
        heatmap(RhoEvolution[:, :, t]', 
            c=:viridis, 
            clims=(min, max), 
            title="Density at step $t", 
            xlabel="X", 
            ylabel="Y", 
            aspect_ratio=1)
    end
    filename = "Animation/DensityEvolution.gif"
    gif(anim, filename, fps=10)
end

end