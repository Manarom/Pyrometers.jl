    """
    signal(p::AbstractPyrometer , Tmeasured::Number)

Returns the signal value which will give the temperature `Tmeasured`
"""
function signal(p::AbstractPyrometer , Tmeasured::Number , ϵ::Number) error("default method") end

signal(p::SingleWavelengthPyrometer , Tmeasured::Number, ϵ::Number)  = ϵ * Planck.ibb(p.λ[] , Tmeasured)
signal(p::SpectralBandPyrometer , Tmeasured::Number, ϵ::Number)= ϵ * Planck.band_power(Tmeasured , λₗ=p.λ[1] , λᵣ=p.λ[2])
ratio_signal(p::TwoWavelengthRatioPyrometer , Tmeasured::Number, ϵ_slope::Number)  = Planck.spectral_ratio(p.λ[1] , p.λ[2] , Tmeasured) * ϵ_slope 
ratio_signal(p::TwoBandsRatioPyrometer , Tmeasured::Number, ϵ_slope::Number) = Planck.spectral_band_ratio(p.λ[1] , p.λ[2] , Tmeasured) * ϵ_slope 

signal(p::TwoBandsRatioPyrometer , Tmeasured::Number) = signal(p , Tmeasured , get_emissivity(p))
signal(p::TwoBandsRatioPyrometer , Tmeasured::Number , ϵ::NTuple{2,D}) where {D <: Number}=begin 
    (b1 , b2) = p.λ[1] , p.λ[2]
     (
        ϵ[1] * Planck.band_power(Tmeasured , λₗ = b1[1] , λᵣ = b1[2]) , 
        ϵ[2] * Planck.band_power(Tmeasured , λₗ = b2[1] , λᵣ = b2[2]) 
    )
end
signal(p::TwoWavelengthRatioPyrometer , Tmeasured::Number) = signal(p , Tmeasured , get_emissivity(p))
signal(p::TwoWavelengthRatioPyrometer , Tmeasured::Number , ϵ::NTuple{2,D}) where D <: Number =begin 
    (l1 , l2) = p.λ[1] , p.λ[2]
     (
        ϵ[1] * Planck.ibb(l1 , Tmeasured ) , 
        ϵ[2] * Planck.ibb(l2 , Tmeasured ) 
    )
end
signal(p::Pyrometer , Tmeasured::Number) = begin
    return signal(p , Tmeasured , _get_epsilon_equivalent(p))
end
signal(p::RatioPyrometer , Tmeasured::Number) = begin
    return signal(p , Tmeasured , get_emissivity(p))
end
# AbstractSpectralQuantity
    """
    signal(p::AbstractPyrometer , Tmeasured::Number, ϵ::AbstractContinuousOrDiscreteQuantity)

Returns the signal value which will give the temperature `Tmeasured`
"""
function signal(p::AbstractPyrometer , Tmeasured::Number, ϵ::AbstractContinuousOrDiscreteQuantity) error(" $(typeof(p)) for $(typeof(ϵ))") end
# single wavelengths and color pyrometers 
signal(p::SingleWavelengthPyrometer , Tmeasured::Number , ϵ::AbstractContinuousOrDiscreteQuantity)  = ϵ(p.λ[] , Tmeasured) * Planck.ibb(p.λ[] , Tmeasured)
function ratio_signal(p::TwoWavelengthRatioPyrometer , Tmeasured::Number, ϵ::AbstractContinuousOrDiscreteQuantity)  
    e1 , e2 = get_single_wavelength_value(p.λ , Tmeasured , ϵ)
    l1 , l2 = p.λ[1] , p.λ[2]
    Planck.spectral_ratio( l1 , l2 , Tmeasured) * e1/e2
end
function signal(p::TwoWavelengthRatioPyrometer , Tmeasured::Number, ϵ::AbstractContinuousOrDiscreteQuantity)  
    e1 , e2 = get_single_wavelength_value(p.λ , Tmeasured , ϵ)
    l1 , l2 = p.λ[1] , p.λ[2]
    (e1 * Planck.ibb(l1 , Tmeasured) , e2 * Planck.ibb(l2 , Tmeasured))
end
#spectral band pyrometer 
signal(p::SpectralBandPyrometer , Tmeasured::Number  , ϵ::AbstractSpectralQuantity) = Planck.planck_weighted(ϵ ,  p.λ[1] , p.λ[2] , Tmeasured)
function signal(p::SpectralBandPyrometer , Tmeasured::Number  , ϵ::AbstractDiscreteQuantity) 
    (l , e) = subrange_view(p.λ[1] , p.λ[2] , ϵ)
    return Planck.planck_weighted(e , l , Tmeasured)
end
#two-bands ratio pyrometer 
ratio_signal(p::TwoBandsRatioPyrometer , Tmeasured::Number, ϵ::AbstractSpectralQuantity) = begin 
    (band1 , band2) = p.λ[1] , p.λ[2]
    return Planck.planck_weighted_ratio(ϵ , band1 , band2 , Tmeasured)
end
signal(p::TwoBandsRatioPyrometer , Tmeasured::Number, ϵ::AbstractSpectralQuantity) = begin 
    (band1 , band2) = p.λ[1] , p.λ[2]
    return (    Planck.planck_weighted(ϵ , band1[1] , band1[2] , Tmeasured) , 
                Planck.planck_weighted(ϵ , band2[1] , band2[2] , Tmeasured) 
            )
end
function ratio_signal(p::TwoBandsRatioPyrometer , Tmeasured::Number  , ϵ::AbstractDiscreteQuantity) 
    ((l1 , l2) , (e1 , e2)) = subrange_view(p.λ[1] , p.λ[2] , ϵ)
    return Planck.planck_weighted_ratio(e1 , l1 , e2 , l2 , Tmeasured)
end
function signal(p::TwoBandsRatioPyrometer , Tmeasured::Number  , ϵ::AbstractDiscreteQuantity) 
    ((l1 , l2) , (e1 , e2)) = subrange_view(p.λ[1] , p.λ[2] , ϵ)
        return (    
                Planck.planck_weighted(e1 , l1 , Tmeasured) , 
                Planck.planck_weighted(e2 , l2 , Tmeasured) 
            )
end