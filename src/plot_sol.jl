import Pkg; Pkg.add("PyPlot")
import Pkg; Pkg.add("Plots")
import PyPlot
using Plots

# Plot one frame
function _plot_diffusion(sol, t::Int, x, y, save; kwargs...)
    # We use a view to avoid doing a copy of the submatrix
    slice = sol[t, :, :]
    plot(x, y, slice; kwargs...)
    if save
        savefig("sol$t.pdf")
    end
end

# Animation
function _plot_diffusion(sol, T::AbstractVector{Int}, x, y, save; kwargs...)
    anim = @animate for t in T
        println("Rendering frame $t/$(length(T))")
        _plot_diffusion(sol, t, x, y, false; kwargs...)

    end
    if save
        gif(anim, "sol.gif", fps=25)
        #mov(anim, "sol.mov", fps=10)
        #mp4(anim, "sol.mp4", fps=10)
    end
end

# sol is an AbstractMatrix{Float64, 3} containing the solution
# If t is an integer, the solution at this time step is saved in "sol$t.png", otherwise, a movie "sol.gif" is saved
function plot_diffusion(sol, t=10:size(sol, 1))
    x = range(0, stop=5, length=size(sol, 2))
    y = range(0, stop=3, length=size(sol, 3))
    zmin = minimum(sol)
    zmax = maximum(sol)
    col = cgrad(:heat, range(zmin, stop=zmax, length=10))
    _plot_diffusion(sol, t, x,y, true, zlims = (zmin, zmax), clims = (zmin, zmax), st = :surface, title = "Temperature profile at t=10s", xlabel = "x   [m]", ylabel = "y   [m]", zlabel = "T   [K]")
    #_plot_diffusion(sol, t, x, y, true, zlims = (-10, 10), clims = (zmin, zmax), st = :surface)
    
end



f = open("sol.txt")
n = map(x -> parse(Int, x), split(readline(f)))
sol = Array{Float64, 3}(undef, n[1], n[2], n[3])
for k in 1:n[1]
    for i in 1:n[3]
        for j in 1:n[2]
            sol[k, j, i] = parse(Float64, readline(f))
        end
    end
end

# Save solution at the 9th timestep in sol9.png
plot_diffusion(sol, 1000)
# Save full solution in sol.gif
#plot_diffusion(sol)
