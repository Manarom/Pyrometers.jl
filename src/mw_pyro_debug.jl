using Revise
using Pkg

Pkg.activate(joinpath(@__DIR__,"..")) 

using Pyrometers

using Optimization , OptimizationOptimJL
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

bbp = Pyrometers.BBPoint(i , l)

Pyrometers.fit_T!(bbp , NelderMead)