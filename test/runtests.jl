using Pyrometers , PlanckFunctions
using ForwardDiff , Zygote , QuadGK , ADTypes
using Test


const pl_fun = PlanckFunctions.ibb

@testset "Pyrometers.jl" begin
    p_vector = Pyrometers.produce_pyrometers() # creating all default pyrometers vector
    p_new = Pyrometers.Pyrometer([1.1; 2.1] , type=:test  , ϵ=0.99) # creating narrow-band pyrometer
    p_new2 = Pyrometers.Pyrometer( 18.0  , type=:test , ϵ=0.99) # creating narrow-band pyrometer
    push!(p_vector,p_new)
    push!(p_vector,p_new2)

    N = length(p_vector)
    Treal = 1235.0 # this is the real temperature of the surface
    Tmeasured = fill(1075.0 , N) #input data, all pyrometers "measure" the same temperature
    Pyrometers.fit_ϵ!(p_vector , 1235.0 , Tmeasured) # fitting the emissivity of each pyrometer

    Tmeasured_calc = Vector{Float64}(undef,N) #
    for (i , p) in enumerate(p_vector)
        @show p
        if Pyrometers.is_spectral_band(p)
            Imeas = PlanckFunctions.band_power(Treal, λₗ=p.λ[1] , λᵣ=p.λ[2] )*p.ϵ[]
        else
            Imeas = PlanckFunctions.ibb(p.λ[1],Treal)*p.ϵ[]
        end
        Tmeasured_calc[i] = p(Imeas)
        @test Tmeasured_calc[i] ≈ Treal  
    end


    l = 0.1:1e-3:10
    a_fun(l) =  0.2 + 0.1*sin(l)
    a =a_fun.(l) 
    l_left = 2.3
    l_right = 4.5
    Treal = 1236.67

    a_analyt(l , t) =  0.5 + sin(l)*1e-4*t + 1e-6*t^2
    d_a_analyt(l , t) = sin(l)*1e-3 + 2e-6*t
    d2_a_analyt(l , t) = 2e-6

    # various formats of emissivity dependent on wavelengh and temperature
    eps_tabular = Pyrometers.TabularQuantity(l , a) # tabular quantity
    eps_isothermal = Pyrometers.IsothermalSpectralQuantity(a_fun) # isothermal quantity 
    eps_analyt = Pyrometers.AnalyticalSpectralQuantity(a_analyt , d_a_analyt , d2_a_analyt) # analytic quantity 
    eps_fdiffed = Pyrometers.GenericDifferentiableSpectralQuantity(a_analyt) # differentiable quantity 

    (l_v , a_v) = Pyrometers.subrange_view(l_left , l_right , eps_tabular)

    # simulated measured signals for spectral band pyrometers 
    imeasure_tuple = ( Pyrometers.Planck.planck_weighted(a_v , l_v  , Treal) , 
                       Pyrometers.Planck.planck_weighted(eps_isothermal , l_left , l_right , Treal), # the same as for the other 
                       first(quadgk(l->eps_analyt(l , Treal)*pl_fun(l , Treal) , l_left , l_right)), 
                       first(quadgk(l->eps_fdiffed(l , Treal)*pl_fun(l , Treal) , l_left , l_right))
                    )

    pyr_band = Pyrometers.SpectralBandPyrometer(l_left , l_right)

    l_single = (l_left + l_right)/2
    pyr_single = Pyrometers.SingleWavelengthPyrometer(l_single)

    i_singl_tup = ( eps_isothermal(l_single , Treal) * pl_fun(l_single , Treal) , 
                    eps_isothermal(l_single , Treal) * pl_fun(l_single , Treal),
                    eps_analyt(l_single , Treal) * pl_fun(l_single , Treal) , 
                    eps_fdiffed(l_single , Treal) * pl_fun(l_single , Treal) 
                    )

    eps_tuple = (eps_tabular , eps_isothermal , eps_analyt , eps_fdiffed)
    for (e , i , i_s) in zip( eps_tuple , imeasure_tuple , i_singl_tup   )
        @show typeof(e)
        @test Pyrometers.measure(pyr_band , i , e) ≈ Treal
        @test Pyrometers.measure(pyr_single , i_s , e) ≈ Treal
    end
    band1 = (2.0 , 3.0)
    band2 = (4.0 , 5.0)
    p_band_ratio = Pyrometers.TwoBandsRatioPyrometer(band1 , band2)
    ((l_v1 , l_v2) ,  (a_v1 , a_v2)) = Pyrometers.subrange_view(band1 , band2 , eps_tabular)

    measured_band_ratio = (Pyrometers.Planck.planck_weighted_ratio(a_v1 , l_v1 , a_v2 , l_v2 , Treal ) , 
                        Pyrometers.Planck.planck_weighted_ratio(eps_isothermal , band1 , band2 , Treal ) ,
                        Pyrometers.Planck.planck_weighted_ratio(eps_analyt , band1 , band2 , Treal ),
                        Pyrometers.Planck.planck_weighted_ratio(eps_fdiffed , band1 , band2 , Treal ))                    
    for (e , i , i_s) in zip( eps_tuple , measured_band_ratio , i_singl_tup   )
        @show typeof(e)
        @test Pyrometers.measure(p_band_ratio , i , e) ≈ Treal
        #@test Pyrometers.measure(pyr_single , i_s , e) ≈ Treal
    end






end

