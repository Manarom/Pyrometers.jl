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
i_int = Pyrometers.integrate(1.0 , 2.0 , Pyrometers.fix_temperature(bb , Ttrue))
Pyrometers.fit_integral(bb ,1.0 , 2.0 , i_int)


bbp = Pyrometers.BBPoint(SVector{50}(i) , SVector{50}(l))
bbp(bb.(l , 1685) )

@benchmark Pyrometers.fit_T!($bbp )
@benchmark $bbp()
@benchmark Pyrometers.fit_T!($bbp , $LBFGS )
@code_warntype Pyrometers.fit_T!(bbp )

bb = PlanckEmitter()
l = range(1,2,50)
Ttrue = 1076.894567 
i = bb.(l , Ttrue)
N = length(l)
noise = 1e-3
mwp = Pyrometers.MWPPoint{N , 3}(i .+ noise* randn(N), SVector{N}(l) , (0.2 , 0.3 , 0.4) , 1273.15)
mwp(i , emissivity_range = (0.5 , 1.0) , e_starting = 0.99)
Pyrometers.clear_cache!(mwp)
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
using StaticArrays
e = Pyrometers.IsothermalSpectralQuantity(l->0.8 + l/10)
bb = Pyrometers.PlanckEmitter()
i = Pyrometers.fix_temperature(e * bb  , 1200.0)
l = range(1,2,50)
N = length(l)
poly_type = Pyrometers.ScaledPolynomials.BernsteinSymPoly{3,Float64}
mwp = Pyrometers.MWPPoint(SVector{N}(i.(l)) , SVector{N}(l) , SVector(0.2 , 0.3 , 0.5) ,   1234.6 , poly_type)
mwp(;emissivity_range = ((0.6 , 0.6 , 0.6) , (0.99 , 0.99 , 0.99)) , temperature_range = (1100.0 , 1300.0) , optimizer = LBFGS)
Pyrometers.emissivity(mwp)
Pyrometers.emissivity_poly(mwp)
Pyrometers.fitting_covariance(mwp)

mwp_pyro = Pyrometers.MultiWavelengthPyrometer{50 , 4}(l ; i_measured = i)
@code_warntype Pyrometers.MultiWavelengthPyrometer{50 , 5}(l ; i_measured = i)
@benchmark Pyrometers.MultiWavelengthPyrometer{50 , 5}($l ; i_measured = $i)

Pyrometers.clear_cache!(mwp_pyro)
mwp_pyro(; emissivity_range = (0.8 , 1.0) , T_starting = 200 , e_starting=0.6)
@benchmark $mwp_pyro()
epoly = Pyrometers.emissivity_poly(mwp_pyro)

poly = poly_type([0.1, 2.0 , 3.4])

mwp_pyro(;optimizer = LBFGS)