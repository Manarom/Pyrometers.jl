using Pyrometers , PlanckFunctions
using ForwardDiff , Zygote , QuadGK , ADTypes
using Test


const pl_fun = PlanckFunctions.ibb
const ratio_fun = PlanckFunctions.spectral_ratio
eps_ratio_fun(e , l1 , l2 , T) = e(l1 , T)/e(l2 , T)

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
        println("Band  and single-wavelength on $(nameof(typeof(e)))")
        @test Pyrometers.measure(pyr_band , i , e) ≈ Treal
        @test Pyrometers.signal(pyr_band , Treal  , e) ≈ i
        @test Pyrometers.measure(pyr_single , i_s , e) ≈ Treal
        @test Pyrometers.signal(pyr_single , Treal  , e) ≈ i_s
    end

    #ratio pyrometer 
    band1 = (2.0 , 3.0)
    band2 = (4.0 , 5.0)
    l_1 , l_2 = sum(band1)/2 , sum(band2)/2 
    p_band_ratio = Pyrometers.TwoBandsRatioPyrometer(band1 , band2)
    pyr_single_ratio = Pyrometers.TwoWavelengthRatioPyrometer(l_1 , l_2)
    ((l_v1 , l_v2) ,  (a_v1 , a_v2)) = Pyrometers.subrange_view(band1 , band2 , eps_tabular)

    measured_band_ratio = (Pyrometers.Planck.planck_weighted_ratio(a_v1 , l_v1 , a_v2 , l_v2 , Treal ) , 
                        Pyrometers.Planck.planck_weighted_ratio(eps_isothermal , band1 , band2 , Treal ) ,
                        Pyrometers.Planck.planck_weighted_ratio(eps_analyt , band1 , band2 , Treal ),
                        Pyrometers.Planck.planck_weighted_ratio(eps_fdiffed , band1 , band2 , Treal ))     
    measured_single_ratio = [ 
        ratio_fun(l_1 , l_2 , Treal) * a_fun(l_1)/a_fun(l_2),
        ratio_fun(l_1 , l_2 , Treal) * eps_ratio_fun(eps_isothermal , l_1 , l_2 , Treal),
        ratio_fun(l_1 , l_2 , Treal) * eps_ratio_fun(eps_analyt , l_1 , l_2 , Treal),
        ratio_fun(l_1 , l_2 , Treal) * eps_ratio_fun(eps_fdiffed , l_1 , l_2 , Treal),
    ]               
    for (e , i , i_s) in zip( eps_tuple , measured_band_ratio , measured_single_ratio   )
        println(" Ratio and band-ratio pyrometers on $(nameof(typeof(e)))")
        @test Pyrometers.measure(p_band_ratio , i , e) ≈ Treal
        @test Pyrometers.signal(p_band_ratio , Treal, e) ≈ i 
        @test Pyrometers.measure(pyr_single_ratio , i_s , e) ≈ Treal
        @test Pyrometers.signal(pyr_single_ratio , Treal , e) ≈ i_s
    end
    # testing the convetr temperature function 
    eps_test = Pyrometers.AnalyticalSpectralQuantity((l , t)-> l + t + t^2  , 
                            (l , t)-> 1 + 2t  , 
                            (l , t)-> 2 )
    
    pyrs = (pyr_band , pyr_single , p_band_ratio , pyr_single_ratio)

    i_test = [ first(quadgk(l->eps_analyt(l , Treal)*pl_fun(l , Treal) , l_left , l_right)),
                    eps_analyt(l_single , Treal) * pl_fun(l_single , Treal) ,
                    Pyrometers.Planck.planck_weighted_ratio(eps_analyt , band1 , band2 , Treal ), 
                    ratio_fun(l_1 , l_2 , Treal) * eps_ratio_fun(eps_analyt , l_1 , l_2 , Treal)
                ]

    t_test  = [ p(i , eps_test)  for (p,i) in zip(pyrs , i_test)]
    for (p,t) in zip(pyrs , t_test)
        @test Pyrometers.convert_temperature(p , t , eps_test , eps_analyt) ≈ Treal
    end

    #testing stray radiation exclution 
    println("Testing stray radiaion correction")
    # planck function integrator within pyrometer's spectral range 
    b_i = Pyrometers.Planck.∫ₗ(Pyrometers.Planck.ibb , 2.0 , 3.0)
    
    T_true = 1200.0  # true temperature of the surface 
    T_src  = 1500.0  # furnace temperature 
    ϵ_src  = 0.85    # furnace emissivity 
    ϵ_surf = 0.7     # surface emissivity 
    i_f(e) = e * b_i(T_true) + (1 - e)*b_i(T_src)
    p = Pyrometers.SpectralBandPyrometer(2.0,3.0, ϵ=ϵ_surf) 

    print("Enclosure geometry ,  fixed emissivity and input signal...")
    geom_enc = Pyrometers.EnclosureGeometry()
    e_eff = Pyrometers.effective_emissivity(geom_enc, ϵ_surf, ϵ_src)
    @test e_eff == ϵ_surf
    T_measured = p(i_f(ϵ_surf))
    T_recovered= Pyrometers.stray_radiation_corrected_temperature(p, T_measured, T_src , ϵ_src , geom_enc)
    @test T_recovered ≈ T_true
    println("ok")

    print("Parallel plates geometry ,  fixed emissivity and input signal...")
    geom = Pyrometers.ViewFactorGeometry(1.0)
    e_eff = Pyrometers.effective_emissivity(geom, ϵ_surf, ϵ_src)
    T_meas = p(i_f(e_eff))
    T_true_par = Pyrometers.stray_radiation_corrected_temperature(p, T_meas, T_src, ϵ_src, geom)
    @test T_true_par ≈ T_true
    println("ok")

   
    print("Fixed surface emissivity incident radiation is provided as a function of wavelength...")
    ϵ_surf = p.ϵ[]
    bb = Pyrometers.PlanckEmitter()
    i_incident = Pyrometers.fix_temperature(bb , 1500.0)
    
    i_full = ϵ_surf * bb + (1 - ϵ_surf) * i_incident 
    i_full_iso = Pyrometers.fix_temperature(i_full , T_true)
    I_total = Pyrometers.integrate(p , i_full_iso)
    T_meas = p(I_total) # measured temperature including stray radiation impact
    # the insident radiation is provided as irradiance 
    T_corrected = Pyrometers.stray_radiation_corrected_temperature(p, T_meas, i_incident) #applying correction 

    @test T_corrected ≈ T_true
    println("ok")

    # ϵ_obj_spec = Pyrometers.AnalyticalSpectralQuantity((λ, t) -> 0.7 - 0.00005 * t , (λ, t) -> - 0.00005   , (λ, t) -> 0.0 )
end

