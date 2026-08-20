using Pyrometers , PlanckFunctions
using ForwardDiff , Zygote , QuadGK , ADTypes
using Test


const pl_fun = PlanckFunctions.ibb

@testset "Pyrometers.jl" begin
    #=p_vector = Pyrometers.produce_pyrometers() # creating all default pyrometers vector
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
        println("Tmeasured_calc[$(i)] = $(Tmeasured_calc[i]); ϵ=$(p.ϵ[])")
        @test Tmeasured_calc[i] ≈ Treal rtol=1e-3 
    end=#


    l = 0.1:1e-3:10
    a_fun(l) =  0.2 + 0.1*sin(l)
    a =a_fun.(l) 
    l_left = 2.3
    l_right = 4.5
    Treal = 1236.67
    # AnalyticalSpectralQuantity
    a_analyt(l , t) =  0.5 + sin(l)*1e-4*t + 1e-6*t^2
    d_a_analyt(l , t) = sin(l)*1e-3 + 2e-6*t
    d2_a_analyt(l , t) = 2e-6
    pyr = Pyrometers.SpectralBandPyrometer(l_left , l_right)
    eps3 = Pyrometers.AnalyticalSpectralQuantity(a_analyt , d_a_analyt , d2_a_analyt )
    imeasure = first(quadgk(l->eps3(l , Treal)*pl_fun(l , Treal) , l_left , l_right))
    eps3(2.0 , 1500)
    @test Pyrometers.measure(pyr , imeasure , eps3) ≈ Treal

    eps4 = Pyrometers.GenericDifferentiableSpectralQuantity(a_analyt , AutoForwardDiff())
    @test Pyrometers.measure(pyr , imeasure , eps4) ≈ Treal

end

