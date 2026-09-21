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

bb = PlanckEmitter()
l = range(1,2,50)
Ttrue = 1076.894567 
i = bb.(l , Ttrue)
N = length(l)
noise = 1e-3
mwp = Pyrometers.MWPPoint(SVector{N}(i .+ 0.0 * randn(N)) , SVector{N}(l) , SVector((0.2 , 0.3 , 0.4 , 1273.15)))
mwp(emissivity_range = (0.1 , 1.0))
Ttrue
@benchmark $mwp(emissivity_range = (0.9 , 1.0))
Pyrometers.fit_T!(mwp , MVector((0.2 , 0.3 , 0.4 , 1273.15)) , nothing ;  emissivity_range = (0.5 , 1.0))      

Pyrometers.evaluate_box_constraints(mwp , (0.8 , 1.0)  )

Pyrometers.robust_lm_search(mwp , MVector(0.2 , 0.9 , 0.4 , 1073.15) , SVector(0.2 , 0.2 , 0.2 , 500.15) , SVector(1.2 , 1.2 , 1.2 , 1473.15))
@benchmark Pyrometers.robust_lm_search($mwp , $(MVector(0.2 , 0.9 , 0.4 , 1073.15)) , $(SVector(0.2 , 0.2 , 0.2 , 500.15)) , $(SVector(1.2 , 1.2 , 1.2 , 1473.15)))

Pyrometers.robust_lm_search(mwp , MVector((0.2 , 0.9 , 0.4 , 1073.15)))
Ttrue

# Создаем базовые константные векторы границ
const low_b = SVector(0.2, 0.2, 0.2, 500.15)
const upp_b = SVector(1.2, 1.2, 1.2, 1473.15)

@benchmark Pyrometers.robust_lm_search(mwp_copy, x_start, $low_b, $upp_b) setup=(
    # Этот блок выполняется перед каждым отдельным тестом:
    mwp_copy = deepcopy($mwp);      # Свежая копия структуры со сброшенными буферами
    # Принудительно портим кэши, чтобы алгоритм гарантированно считал заново
    mwp_copy.x_jac_vec .= Inf;
    mwp_copy.x_hess_vec .= Inf;
    mwp_copy.x_hess_approx .= Inf;
    
    x_start = MVector(0.2, 0.9, 0.4, 1073.15); # Свежая стартовая точка
)

e = Pyrometers.IsothermalSpectralQuantity(l->0.8 + l/10)
bb = PlanckEmitter()
i = Pyrometers.fix_temperature(e * bb  , 1200.0)
N = length(l)
mwp = Pyrometers.MWPPoint(SVector{N}(i.(l)) , SVector{N}(l) , SVector(0.2 , 0.3 , 0.5 , 1234.6))
mwp(;emissivity_range = ((0.6 , 0.6 , 0.7) , (0.8 , 0.8 , 0.8)) , temperature_range = (1100.0 , 1300.0))
Pyrometers.emissivity(mwp)
