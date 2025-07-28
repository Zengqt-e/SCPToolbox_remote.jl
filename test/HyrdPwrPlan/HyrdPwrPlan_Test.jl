import Pkg
Pkg.activate(".")

using SCPToolbox
using PyPlot, Colors, LinearAlgebra

using ECOS

pbm = TrajectoryProblem();

## set the system dynamics and its derivative
n, m, d = 1, 2, 1;  # set the system dimensions
problem_set_dims!(pbm, n, m, d);

t_f = 0.1;  # set the regular planning period 100ms
Cap = 450e-6;   #450uF capacitor
LmdG = 0.06;    # 0.06Wb
weG = 3e3/60*4*2*pi; #3krpm, 1256.6
wmM = 1e3/60*2*pi;   #1krpm, 104.72

f(t, x, u, p) = begin       # c is constant/parameters in system
    udc, = x;
    iqG, tqM = u;
    return [-1/450e-6*(1.5*iqG*1256.6*0.06 + tqM*104.72)/udc*t_f];
end

A(t, x, u, p) = begin    # A*δx
    udc, = x;
    iqG, tqM = u;
   return [-1/450e-6*(1.5*iqG*1256.6*0.06 + tqM*104.72)*(-1)*udc^(-2)*t_f;;];
end

B(t, x, u, p) =begin     # B*δu = B*δ(iqG;tqM)
    udc, = x;
    iqG, tqM = u;
    return  -1/450e-6*[1.5*1256.6*0.06  104.72]*udc^(-1)*t_f;
end

F(t, x, u, p) =begin
    return zeros(1,1);
end

wrap(func) = (t, k, x, u, p, pbm) -> func(t, x, u, p);
problem_set_dynamics!(pbm, wrap(f), wrap(A), wrap(B), wrap(F));

## set the boundary conditions 
# from 300V to 400V in 100ms
x_0 = [300];
x_f = [400];

g_ic(x, p) = x - x_0;
g_tc(x, p) = x - x_f;
H_0(x, p) = I(1);
H_f(x, p) = I(1);       # 除了状态，还需要一阶导的终止条件

wrap(func) = (x, p, pbm) -> func(x, p);
problem_set_bc!(pbm, :ic, wrap(g_ic), wrap(H_0));
problem_set_bc!(pbm, :tc, wrap(g_tc), wrap(H_f));

## set the desired algorithm
alg = :ptr;
N, Nsub = 20, 25

## set objective function
tqM_req = straightline_interpolate([0], [100], N);;   # 50Nm output is desired
#iqG_Fwd = -tqM_req*104.72/(1.5*1256.6*0.06); 
iqG_Fwd = -tqM_req;
R(x, u ,p) = [(u[2]-50)]'*[(u[2]-50)];
wrap(func) = (t, k, x, u, p, pbm) -> func(x, u, p);
problem_set_running_cost!(pbm, alg, wrap(R));

## set the initial guess
state_guess(N)=straightline_interpolate(x_0, x_f, N);
input_guess(N)= straightline_interpolate([iqG_Fwd[1], tqM_req[1]], [iqG_Fwd[N], tqM_req[N]], N);
problem_set_guess!(pbm, (N, pbm) -> begin
    x = state_guess(N)
    u = input_guess(N)
    p = zeros(1)
    return x, u, p
end)

# PTR Parameters
N, Nsub = 20, 25
iter_max = 50
disc_method = FOH
wvc, wtr = 1e3, 1e0
feas_tol = 1e-2
ε_abs, ε_rel = 1e-3, 1e-2
q_tr = Inf
q_exit = Inf
solver, solver_options = ECOS, Dict("verbose"=>0)

pars = PTR.Parameters(N, Nsub, iter_max, disc_method, wvc, wtr, ε_abs,
                      ε_rel, feas_tol, q_tr, q_exit, solver, solver_options);

if alg == :ptr
    ptr_pbm = PTR.create(pars, pbm)
    sol, history = PTR.solve(ptr_pbm)
end;

sol_PwrPlan = sol;