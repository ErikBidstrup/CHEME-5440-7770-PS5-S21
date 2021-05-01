# Example for constructing a streamplot in Julia for phase portraits
# CHEME5440/7770 - Spring 2021
# The example makes use of Makie.jl and AbstractPlotting.jl for plotting
# To install the packages, use the following commands inside the Julia REPL:
# using Pkg
# Pkg.add("Makie")
# Pkg.add("AbstractPlotting")

using Makie
using AbstractPlotting.MakieLayout
using AbstractPlotting


AbstractPlotting.inline!(true)


# Model for precise adaptation
# D1: active delta receptor in cell 1
# D2: active delta receptor in cell 2
function precise_adapt(D1, D2)

    f_D1 = (D1.^2) ./ (0.1+D1.^2)
    f_D2 = (D2.^2) ./ (0.1+D2.^2)

    u = 1 ./ (1+10*(f_D2).^2) - D1     #dD1/dtau*
    v = 1 ./ (1+10*(f_D1).^2) - D2     #dD2/dtau*
    
    return Point(u,v)
end

# Construct the streamplot
plt1, layout = layoutscene(resolution = (1000,800))

ax = layout[1, 1] = Axis(plt1, xlabel = "D1", ylabel = "D2",
    title = "Streamplot of D1 vs. D2")

streamplot!(ax,precise_adapt, 0..1, 0..1, colormap = :plasma, 
    gridsize= (32,32), arrow_size = 0.025)

# Plotting nullclines for dD1/dtau* = 0 and dD2/dtau* = 0
D1 = LinRange(0,1,100)
f_D1 = (D1.^2) ./ (0.1 .+ D1.^2)
D2 = 1 ./ (1 .+ 10 .* (f_D1).^2)
lines!(ax,D1,D2, linestyle = "-", label  = "dD1/dtau* = 0", linewidth = 5)
lines!(ax,D2,D1, color = :red, linestyle = "-", label  = "dD2/dtau* = 0", linewidth = 5)
layout[1, 2] = Legend(plt1, ax, "Nullclines", framevisible = false)

# Display the plot
display(plt1)

# Save the plot
save("odeField_1b.png")