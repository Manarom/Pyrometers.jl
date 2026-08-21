using BenchmarkTools
using Revise
using QuadGK
using Test
using ForwardDiff
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
@benchmark Pyrometers.fit_ϵ($(p_vector[1]), $Tmeasured[1] , $Treal)

# testing discrete data integration and coninuous function wrapper
#begin 
    l = 0.1:1e-3:10
    a_fun(l) =  0.2 + 0.1*sin(l)
    a =a_fun.(l) 
    l_left = 2.3
    l_right = 4.5
    Treal = 1236.67

    a_analyt(l , t) =  0.5 + sin(l)*1e-4*t + 1e-6*t^2
    d_a_analyt(l , t) = sin(l)*1e-3 + 2e-6*t
    d2_a_analyt(l , t) = 2e-6


    eps_tabular = Pyrometers.TabularQuantity(l , a)
    eps_isothermal = Pyrometers.IsothermalSpectralQuantity(a_fun)
    eps_analyt = Pyrometers.AnalyticalSpectralQuantity(a_analyt , d_a_analyt , d2_a_analyt)
    
    (l_v , a_v) = Pyrometers.subrange_view(l_left , l_right , eps_tabular)
    
    pl_fun(l ,t) = Pyrometers.Planck.ibb(l , t)


    imeasure_tuple = ( Pyrometers.Planck.planck_weighted(a_v , l_v  , Treal) , 
                       Pyrometers.Planck.planck_weighted(eps_isothermal , l_left , l_right , Treal),
                       first(quadgk(l->eps_analyt(l , Treal)*pl_fun(l , Treal) , l_left , l_right)) )

    pyr_band = Pyrometers.SpectralBandPyrometer(l_left , l_right)

    l_single = (l_left + l_right)/2
    pyr_single = Pyrometers.SingleWavelengthPyrometer(l_single)
    i_singl_tup = (eps_isothermal(l_single , Treal) * pl_fun(l_single , Treal) , 
               eps_isothermal(l_single , Treal) * pl_fun(l_single , Treal),
               eps_analyt(l_single , Treal) * pl_fun(l_single , Treal))
    ctx_i = nothing 
    ctx_is = nothing
    eps_tuple = (eps_tabular , eps_isothermal , eps_analyt)
    for (e , i , i_s) in zip( eps_tuple , imeasure_tuple , i_singl_tup   )
        @show typeof(e)
        @test Pyrometers.measure(pyr_band , i , e) ≈ Treal
        @test Pyrometers.measure(pyr_single , i_s , e) ≈ Treal
    end
    eps_fdiffed = Pyrometers.GenericDifferentiableSpectralQuantity(a_analyt)
    ctx = Pyrometers.GenericSpectralQuantityContext(pyr_band , eps_fdiffed , imeasure_tuple[3])
    ctx(Treal)



band1 = (2.0 , 3.0)
band2 = (4.0 , 5.0)
p_band_ratio = Pyrometers.TwoBandsRatioPyrometer(band1 , band2)
((l_v1 , l_v2) ,  (a_v1 , a_v2)) = Pyrometers.subrange_view(band1 , band2 , eps_tabular)

measured_band_ratio = (Pyrometers.Planck.planck_weighted_ratio(a_v1 , l_v1 , a_v2 , l_v2 , Treal ) , 
                       Pyrometers.Planck.planck_weighted_ratio(eps_isothermal , band1 , band2 , Treal ) ,
                       Pyrometers.Planck.planck_weighted_ratio(eps_analyt , band1 , band2 , Treal ),
                       Pyrometers.Planck.planck_weighted_ratio(eps_fdiffed , band1 , band2 , Treal ))                    
eps_tuple[1]
Pyrometers.measure(p_band_ratio , measured_band_ratio[1] , eps_tuple[1])
    for (e , i , i_s) in zip( eps_tuple , measured_band_ratio , i_singl_tup   )
        @show typeof(e)
        @test Pyrometers.measure(p_band_ratio , i , e) ≈ Treal
        #@test Pyrometers.measure(pyr_single , i_s , e) ≈ Treal
    end







ctx = Pyrometers.GenericSpectralQuantityContext(p_band_ratio , eps_tabular , measured_ratio)
ctx(1236.67)
Pyrometers.measure(p_band_ratio , measured_ratio , eps_tabular)


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
ctx = Pyrometers.TabularQuantityContext(p1 , eps , i_meas)
ctx(1200)
Pyrometers.measure(p1  , i_meas, eps)
Pyrometers.measure(p1  , i_meas , eps2)
typeof(eps)

p2 = Pyrometers.TwoWavelengthRatioPyrometer(2.0 ,3.0)
