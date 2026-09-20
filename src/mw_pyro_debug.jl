using Revise

using Pkg

Pkg.activate(joinpath(@__DIR__,"..")) 

using Pyrometers
using StaticArrays
#using Optimization , OptimizationOptimJL
using Optim
#using JSOSolvers
using BenchmarkTools
using QuadGK
using Test
using ForwardDiff

using ADTypes
using PlanckFunctions , DataInterpolations
using StaticArrays
bb = PlanckEmitter()
l = range(1,2,50)
Ttrue = 1076.894567 
i = bb.(l , Ttrue)

bbp = Pyrometers.BBPoint(SVector{50}(i) , SVector{50}(l))
bbp(bb.(l , 1685) )

@benchmark Pyrometers._solve_problem($bbp , $[1237.8] , $[20.0] , $[3000.0] , LBFGS)

@benchmark $bbp(;optimizer=LBFGS)
@benchmark Pyrometers.fit_T!($bbp )
@benchmark $bbp()
@benchmark Pyrometers.fit_T!($bbp , $LBFGS )
@code_warntype Pyrometers.fit_T!(bbp )

optim_fun = OptimizationFunction(Pyrometers.disc , grad = Pyrometers.grad! , hess = Pyrometers.hess!)
optim_fun([1274.0] , bbp)
prob = OptimizationProblem(optim_fun , [1000.0] , bbp)
solve(prob , LBFGS())
starting_vector = [1000.0]
probl = OptimizationProblem(
                    optim_fun, 
                    starting_vector,
                    bbp, 
                )
solve(probl , LBFGS())
@benchmark solve($probl , $(Brent()))

N = length(l)
noise = 1e-3
mwp = Pyrometers.MWPPoint(SVector{N}(i .+ 0.0 * randn(N)) , SVector{N}(l) , SVector((0.2 , 0.3 , 0.4 , 1273.15)))
mwp(emissivity_range = (0.9 , 1.0))
Ttrue
@benchmark $mwp()
Pyrometers.fit_T!(mwp , MVector((0.2 , 0.3 , 0.4 , 1273.15)) , nothing ;  emissivity_range = (0.5 , 1.0))      

Pyrometers.evaluate_box_constraints(mwp , (0.2 , 0.3)  , (1200.0 , 1300.0) )