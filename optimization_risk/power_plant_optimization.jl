using DataFrames, LinearAlgebra, Statistics, Distributions
using JuMP, Ipopt, HiGHS, GLPK

N = 3                        # Number of plants
F = [3000.0, 1200.0, 2500.0] # Fixed cost $
V = [5.0,7.0,3.0]            # Variable cost $/MWh   
C = [7000, 800.0, 3000.0]    # Max production capacity MWh
D = 8000.0                   # Total domand MWh
B = 7.0                      # $/MWh from external source
S = 10.0                     # Selling price $/MWh
Mi = 1000                    # Maximum import 1000 MWh

model = Model(GLPK.Optimizer)
set_silent(model)
@variable(model, x[1:N], Bin)
@variable(model, q[1:N] >= 0)
@variable(model, e >= 0)
@constraint(model, sum(q[i] for i in 1:N) + e == D)
@constraint(model, e <= Mi)
for i in 1:N
    @constraint(model, q[i] <= C[i] * x[i])
end
@objective(model, Max, 
    S * (sum(q[i] for i in 1:N) + e) - 
    (sum(F[i] * x[i] for i in 1:N) + sum(V[i] * q[i] for i in 1:N) + B * e)
)
optimize!(model)
value.(x)
value.(e)
value.(q)
