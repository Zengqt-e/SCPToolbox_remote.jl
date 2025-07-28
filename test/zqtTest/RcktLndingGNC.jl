import Pkg
import Pkg
Pkg.activate("..")

## these lines are required for local installations
# Pkg.develop(path="../../SCPToolbox.jl/")
# Pkg.precompile()

using SCPToolbox
using PyCall, PyPlot, Colors, LinearAlgebra

# Import the different possible low-level convex solvers
#using ECOS
# using Ipopt
# using SCS
# using OSQP

pbm = TrajetoryProblem();

# Environment parameters
g = 1.625;       #[m/s^2]
g_E = 9.807;     #[m/s^2]
# Mechanical parameters
m_wet = 25e3;   # [kg] Initial mass
L = 0.5;        # [m] Thrust lever arm
J = 100e3;      # [kg*m^2] Moment of inertia
# Propulsion parameters
Isp = 370;      # [s] Specific Impulse 
T_min = 20e3;   # [N] Minimum Thrust
T_max = 80e3;   # [N] Maximum Thrust
α = 1/(Isp*g_E);    # [kg/s/N]
δ




