# Example for use of DifferentialEquations.jl
# CHEME5440/7770 - Spring 2021
# To install the package, use the following commands inside the Julia REPL:
# using Pkg
# Pkg.add("DifferentialEquations")
# For this example, you will also need Plots; To add:
# Pkg.add("Plots")

# ------------------------------------------------------------ #
# Lorenz system example
# dX/dt = σ(Y-X)
# dY/dt = X(ρ - Z) - Y
# dZ/dt = XY - βZ
# Based on DifferentialEquations.jl tutorial

using DifferentialEquations     # Include DifferentialEquations.jl
using Plots                     # Include Plots.jl for plotting
gr(show = true)  # Use the gr backend for plotting and show plots


# Model parameters:             Note: 0 indicates to the compiler that ρ is a floating-point number (i.e. not an integer) 
lig = 10.0 # nM
k_r = 0.1 # s^-1
E0₀ = 3500.0 # nM
v_rmax = k_r*E0₀ # nM s^-1
a_b = 1.0 # s^-1 nM^-1
a_bp = 0.1 # s^-1 nM^-1 .
d_b = 1.0 # s^-1
d_bp = 0.01 # s^-1
k_b = 0.0 # s^-1
k_bp = 1.0 # s^-1
k_1pos = (1.0/1000.0) #s^-1 nM^-1
k_1neg = 1.0 #s^-1
a_1pos = (1) / (1.0+lig) # s^-1
a_1neg = (lig) / (1.0+lig) # s^-1
beta_1 = (2.5*lig) / (1.0+lig) # s^-1


# -------------------------- Lorenz Model -------------------------------------
# du: Differential equations, [dX/dt, dY/dt, dZ/dt]
# u: Time-dependent variables, [X, Y, Z]; u[1] = xX; u[2] = Y; u[3] = Z
# p: Additional model parameters (none in this example)
# t: time 
# Note "!" point after function name is a Julia convention that indicates 
# that the function will modify values in one or more of the function arguments.  In this case,
# lorenz! will modify values in the input vector, du
function lorenz!(du,u,p,t)
    du[1] = -v_rmax + k_bp*u[7] + k_b*u[6]
    du[2] = v_rmax + a_1neg*u[3] - a_1pos*u[2] + beta_1*(u[6]+u[7])
    du[3] = a_1pos*u[2] - a_1neg*u[3] - a_bp*u[3]*u[5] - a_b*u[3]*u[4] + d_bp*u[7] + d_b*u[6]
    du[4] = -k_1pos*u[4]*u[3] + k_1neg*u[5] + beta_1*u[6] + k_b*u[6] - a_b*u[4]*u[3] + d_b*u[6]
    du[5] = k_1pos*u[4]*u[3] - k_1neg*u[5] + beta_1*u[7] + k_bp*u[7] - a_bp*u[5]*u[3] + d_bp*u[7]
    du[6] = a_b*u[3]*u[4] - d_b*u[6] - k_b*u[6] - beta_1*u[6]
    du[7] = a_bp*u[3]*u[5] - d_bp*u[7] - k_bp*u[7] - beta_1*u[7]
end

# ------------- SOLVE THE MODEL WITH DIFFERENTIALEQUATIONS.jl -------------------
E0₀ = 5507.0              # initial value of E0
E1₀ = 350.0                # initial value of E1
E1_star₀ = 2146.0          # initial value of E1*
B₀ = 0.768              # initial value of B
Bp₀ = 1.65                # initial value of Bp
E1_B₀ = 1648.0              # initial value of {E1*B}
E1_Bp₀ = 350.0            # initial value of {E1*Bp}

u0 = [E0₀; E1₀; E1_star₀; B₀; Bp₀; E1_B₀; E1_Bp₀]       # initial state vector
tspan = (0.0,5000)                     #time interval (start time, end time)
prob = ODEProblem(lorenz!,u0,tspan)     #Create an ODE problem for the Lorenz fxn
sol = solve(prob)                       #Solve the system

# ------------- MAKE SOME PLOTS OF THE RESULTS WITH PLOTS.jl ---------------------
#Plot the results; X, Y, and Z vs time
plt1 = plot(sol, xaxis="time", yaxis = "Conc (nM)", label=["E0" "E1" "E1*" "B" "B_p" "E*B" "E*B_p"])
display(plt1)

#Plot the results; the vars=(0,3) argument specifies to plot Z (column 3 of sol)
#vs t
plt2 = plot(sol[3,:]/E1_star₀, xaxis="time", yaxis = "A/A_st", legend = false)
display(plt2)

#Save the plot as PNG files
savefig(plt2, "./plot.png")
savefig(plt1, "./plot1.png") 