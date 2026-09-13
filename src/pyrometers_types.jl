
"""
    SingleWavelengthPyrometer{N, T, DT} <: AbstractPyrometer{N, T}

Type single -wavelengh pyrometer system. 

# Fields
- `type::Symbol`: A unique identifier or descriptive name for the pyrometer model.
- `λ::SVector{1, DT}`: working wavelength ,`μm`
- `ϵ::Base.RefValue{T}`: In-place mutable emissivity value.

# Constructors
    SingleWavelengthPyrometer(λ::Number; type::Symbol=:def, ϵ::Number=1.0)
# Examples
```julia
# Single-wavelength pyrometer (0.85 μm)
p = SingleWavelengthPyrometer(0.85, ϵ=0.9)
```
"""
const SingleWavelengthPyrometer = Pyrometer{1}
"""
    SpectralBandPyrometer

Type alias for a brightness pyrometer over **a defined spectral band**. 
Evaluates signals by integrating the Planck distribution function over the defined spectral range.
    
# Fields
- `type::Symbol`: A unique identifier or descriptive name for the pyrometer model.
- `λ::SVector{2, DT}`: Spectral range 
- `ϵ::Base.RefValue{T}`: In-place mutable emissivity value.

# Constructors
    SpectralBandPyrometer(λ1::Number , λ2::Number ; type::Symbol=:def, ϵ::Number=1.0)
    SpectralBandPyrometer( λ1::NTuple{2 , <:Number} ; type::Symbol=:def, ϵ::Number=1.0)
# Examples
```julia
# Spectral band pyrometer for 2.0 - 4.5 μm spectral band 
p_band = SpectralBandPyrometer(2.0 , 4.5 , ϵ=0.9)
```
"""
const SpectralBandPyrometer = Pyrometer{2}
#Constructors

SingleWavelengthPyrometer(λ::Number; type::Symbol=:def, ϵ::Number=1.0) = Pyrometer(λ; type=type, ϵ=ϵ)

SpectralBandPyrometer(λₗ::Number, λᵣ::Number; type::Symbol=:def, ϵ::Number=1.0) = Pyrometer((λₗ, λᵣ); type=type, ϵ=ϵ)

SpectralBandPyrometer(λ::NTuple{2}; type::Symbol=:def, ϵ::Number=1.0) = Pyrometer(λ ; type=type, ϵ=ϵ)
"""
    TwoWavelengthRatioPyrometer

Two-color spectral ratio pyrometer, working on two single wavelengths 
    
# Fields
- `type::Symbol`: A unique identifier or descriptive name for the pyrometer model.
- `λ1::Number`: First wl
- `λ1::Number`: Second wl
- `ϵ1::Base.RefValue{T}`: First wavelength emissivity.
- `ϵ2::Base.RefValue{T}`: Second wavelength emissivity.

# Constructors
    TwoWavelengthRatioPyrometer(λ1::Number , λ2::Number ; type::Symbol=:def, ϵ1::Number=1.0 , ϵ2::Number=1.0 )
    TwoWavelengthRatioPyrometer( λ::NTuple{2 , <:Number} ; type::Symbol=:def, ϵ1::Number=1.0 , ϵ2::Number=1.0)
# Examples
```julia
# Two-color pyrometer for 2.0 and  4.5 μm wavelength
p_band = TwoWavelengthRatioPyrometer(2.0 , 4.5 , ϵ1=0.9 , e2 = 0.5)
```
"""
const TwoWavelengthRatioPyrometer = RatioPyrometer{2 , T , DT} where {T , DT <: Number}
"""
    TwoBandsRatioPyrometer

Two bands spectral ratio pyrometer, working on two wide spectral bands 
    
# Fields
- `type::Symbol`: A unique identifier or descriptive name for the pyrometer model.
- `λ1::NTuple{2, <:Number}`: First band
- `λ1::NTuple{2, <:Number}`: Second band
- `ϵ1::Base.RefValue{T}`: First band emissivity.
- `ϵ2::Base.RefValue{T}`: Second band emissivity.

# Constructors
    TwoBandsRatioPyrometer(λ1::NTuple{2 , <:Number} , λ2::NTuple{2 , <:Number} ; type::Symbol=:def, ϵ::Number=1.0)
# Examples
```julia
# Spectral band pyrometer for 2.0 - 4.5 μm spectral band 
p_band = TwoBandsRatioPyrometer((2.0 , 4.5) , (7.0 , 9.1) , ϵ1 = 0.4 , ϵ2 = 0.93)
```
"""
const TwoBandsRatioPyrometer = RatioPyrometer{2 , T , DT} where {T , DT <: Tuple}

# Constructor for TwoWavelengthRatioPyrometer (discrete wavelengths)
TwoWavelengthRatioPyrometer(λ1::Number, λ2::Number; type::Symbol=:def, ϵ1::Number=1.0, ϵ2::Number=1.0) = RatioPyrometer((λ1, λ2); type=type, ϵ1=ϵ1, ϵ2=ϵ2)

TwoWavelengthRatioPyrometer(λ::NTuple{2 , D}; type::Symbol=:def, ϵ1::Number=1.0, ϵ2::Number=1.0) where D<: Number = RatioPyrometer(λ; type=type, ϵ1=ϵ1, ϵ2=ϵ2)

# Constructor for TwoBandsRatioPyrometer (integrated spectral bands)
TwoBandsRatioPyrometer(band1::NTuple{2, <:Number}, band2::NTuple{2, <:Number}; 
                               type::Symbol=:def, ϵ1::Number=1.0, 
                               ϵ2::Number=1.0) = RatioPyrometer((band1, band2); type=type, ϵ1=ϵ1, ϵ2=ϵ2)

    """
    wlength(::Pyrometer{N}) where N

Returns the number of wavelengths 
"""
wlength(::Pyrometer{N}) where N = N
"""
    is_spectral_band(p::AbstractPyrometer)

True if `p` wprks in spectral band (not at single wavelength(s))
"""
is_spectral_band(::AbstractPyrometer)  = false
is_spectral_band(::SpectralBandPyrometer)  = true
is_spectral_band(::TwoBandsRatioPyrometer) = true
"""
    is_single_wavelength(::Pyrometer{N}) where N

True if the Pyrometer is a single wavelength or two-color
"""
is_single_wavelength(p::AbstractPyrometer) = ~is_spectral_band(p)

is_spectral_ratio(::AbstractPyrometer) = false 
is_spectral_ratio(::RatioPyrometer) = true

    """
    _get_epsilon_equivalent(p::AbstractPyrometer)

Returns emissivity equivalent to be used in computations, for partial radiation pyrometers 
`typeof(p) <: Pyrometer` returns emissivity value , for  `typeof(p) <: RatioPyrometer` gives
`e-slope` value (two-bands equivalent emissivities ratio )
"""
_get_epsilon_equivalent(p::AbstractPyrometer)= error("not implemented")

_get_epsilon_equivalent(p::Pyrometer)  = p.ϵ[]

_get_epsilon_equivalent(p::RatioPyrometer)  = e_slope(p)

get_emissivity(p::AbstractPyrometer) = error("not implemented")
get_emissivity(p::Pyrometer) = p.ϵ[]
get_emissivity(p::RatioPyrometer) = (p.ϵ1[] , p.ϵ2[])

extract_pyrometer_inds(p::SpectralBandPyrometer , λ::AbstractVector ) =    extract_subrange_inds(p.λ[1] , p.λ[2] , λ)
extract_subrange_inds(l1 , l2 , λ) = (searchsortedfirst(λ , l1 ) , searchsortedlast( λ , l2))