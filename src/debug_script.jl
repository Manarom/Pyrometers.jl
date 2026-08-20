using BenchmarkTools
using Revise
using QuadGK
includet("Pyrometers.jl")
Pyrometers.switch_the_type(8.5)

# default pyrometers testing
p_vector = Pyrometers.produce_pyrometers() # creating all default pyrometers vector
p_new = Pyrometers.Pyrometer([1.1; 2.1] , type=:test,ϵ=0.99) # creating narrow-band pyrometer
p_new2 = Pyrometers.Pyrometer(18.0 , type=:test , ϵ=0.99) # creating narrow-band pyrometer
push!(p_vector,p_new)
push!(p_vector,p_new2)
N = length(p_vector)
Treal = 1435.0 # this is the real temperature of the surface
Tmeasured = fill(1375.0,N) #input data, all pyrometers "measure" the same temperature
Pyrometers.fit_ϵ.(p_vector, Tmeasured , Treal) # fitting the emissivity of each pyrometer
@benchmark Pyrometers.fit_ϵ!($p_vector,1335.0,$Tmeasured) 
@benchmark Pyrometers.fit_ϵ($p_vector[1], $Tmeasured[1] , Treal)

# testing discrete data integration and coninuous function wrapper
l = 0.1:1e-3:10
a_fun(l) =  0.2 + 0.1*sin(l)
a =a_fun.(l) 
l_left = 2.3
l_right = 4.5
Treal = 1236.67
f = @. (l <= l_right) & (l >= l_left)
a_v = @view a[f]
l_v = @view l[f]
i_meas = first(Pyrometers.Planck.planck_weighted(a_v , l_v  , Treal))

# TabularQuantity
pyr = Pyrometers.SpectralBandPyrometer(l_left , l_right)
eps = Pyrometers.TabularQuantity(l , a)
ctx = Pyrometers.TabularQuantityContext(pyr.λ , eps , i_meas)
ctx(1236.67)
Pyrometers.measure(pyr , i_meas , eps)
t = @btime Pyrometers.measure($pyr , i_meas , $eps)
@benchmark Pyrometers.measure($pyr , i_meas , $eps)

# IsothermalSpectralQuantity
eps2 = Pyrometers.IsothermalSpectralQuantity(a_fun)
eps2(2.0 , 1200)
pl_fun = Pyrometers.Planck.ibb
imeasure = quadgk(l->eps2(l , Treal)*pl_fun(l , Treal) , l_left , l_right)
Pyrometers.measure(pyr , imeasure[1] , eps2)
@benchmark Pyrometers.measure($pyr , $imeasure[1] , $eps2)

# AnalyticalSpectralQuantity
a_analyt(l , t) =  0.5 + sin(l)*1e-4*t + 1e-6*t^2
d_a_analyt(l , t) = sin(l)*1e-3 + 2e-6*t
d2_a_analyt(l , t) = 2e-6
eps3 = Pyrometers.AnalyticalSpectralQuantity(a_analyt , d_a_analyt , d2_a_analyt )
imeasure = first(quadgk(l->eps3(l , Treal)*pl_fun(l , Treal) , l_left , l_right))
eps3(2.0 , 1500)
Pyrometers.measure(pyr , imeasure , eps3)

p2 = Pyrometers.TwoBandsRatioPyrometer((2.0 , 3.0) , (4.0 , 5.0))
f1 = @. (l <= 3.0) & (l >= 2.0 )
f2 = @. (l <= 5.0) & (l >= 4.0)
a_v1 = @view a[f1]
l_v1 = @view l[f1]
a_v2 = @view a[f2]
l_v2 = @view l[f2]

i_meas1 = Pyrometers.Planck.planck_weighted(a_v1 , l_v1  , Treal )
i_meas2 = Pyrometers.Planck.planck_weighted(a_v2 , l_v2  , Treal )

i_meas_ratio = i_meas1/ i_meas2
ctx = Pyrometers.TabularQuantityContext(p2.λ[1] , p2.λ[2] , eps , i_meas_ratio)
ctx(1236.67)
Pyrometers.measure(p2 , i_meas_ratio , eps)
Pyrometers.measure(p2 , i_meas_ratio , eps2)

@benchmark Pyrometers.measure($p2 , $i_meas_ratio , $eps)
@benchmark Pyrometers.measure($p2 , $i_meas_ratio , $eps2)

cntx2 = Pyrometers.GenericSpectralQuantityContext(p2.λ , eps3 , i_meas_ratio)
@benchmark $cntx2( 1300)
cntx2(1236.67)
eps4 = Pyrometers.GenericDifferentiableSpectralQuantity(a_analyt)
eps4(2.0 , 1900.0)


p1 = Pyrometers.SingleWavelengthPyrometer(2.0)
Pyrometers.measure(p1  , i_meas, eps)
Pyrometers.measure(p1  , i_meas , eps2)
typeof(eps)

p2 = Pyrometers.TwoWavelengthRatioPyrometer(2.0 ,3.0)
