

"""
    Pyrometer type implements the following methods:

    [``]

"""
    abstract type AbstractPyrometer{N , T} end

"""
    RatioPyrometer{N, T, DT} <: AbstractPyrometer{N, T}

Type representing a ratio (two-color) pyrometer system. 

Supports both two-wavelength monochromatic systems and two-band ratio systems.

# Fields
- `type::Symbol`: A unique identifier or descriptive name for the pyrometer model.
- `λ::NTuple{2, DT}`: The spectral parameters. For monochromatic systems, `DT <: Number` (a pair of wavelengths). For band-ratio systems, `DT <: Tuple` (a pair of wavelength ranges `(λₗ, λᵣ)`).
- `ϵ1::Base.RefValue{T}`: In-place mutable emissivity value for the first channel/band.
- `ϵ2::Base.RefValue{T}`: In-place mutable emissivity value for the second channel/band.

# Constructors
    RatioPyrometer(λ::Union{NTuple{2, T}, NTuple{2, NTuple{2, T}}}; type::Symbol=:def, ϵ1::Number=1.0, ϵ2::Number=1.0)

Construct a `RatioPyrometer` object. The type parameter is automatically inferred as either monochromatic or band-based depending on the structure of `λ`.

# Examples
```julia
# Two-color monochromatic pyrometer (0.85 μm and 1.05 μm)
p_mono = RatioPyrometer((0.85, 1.05), ϵ1=0.9, ϵ2=0.88)

# Two-band ratio pyrometer (bands: 1.5–1.8 μm and 2.0–2.4 μm)
p_band = RatioPyrometer(((1.5, 1.8), (2.0, 2.4)), ϵ1=0.5, ϵ2=0.5)
```
"""
    struct RatioPyrometer{N , T , DT} <: AbstractPyrometer{N , T}
        type::Symbol 
        λ::NTuple{2 , DT}
        ϵ1::Base.RefValue{T}
        ϵ2::Base.RefValue{T}
        RatioPyrometer(λ::Union{NTuple{2 , T} , NTuple{2 , NTuple{2,T}}}; 
                            type::Symbol=:def ,  
                            ϵ1::Number=1.0 , ϵ2::Number = 1.0) where T <: Number= begin
            new{2 , T , eltype(λ)}(type , λ , Ref(ϵ1) ,  Ref(ϵ2))
        end
    end
    e_slope(p::RatioPyrometer) = p.ϵ1[]/p.ϵ2[]
"""
    Pyrometer{N, T} <: AbstractPyrometer{N, T}

Type representing a radiation (brightness) pyrometer system.

Supports both single-wavelength (monochromatic) and narrow-band (integrated spectral range) instruments.

# Fields
- `type::Symbol`: A unique identifier or descriptive name for the pyrometer model.
- `λ::SVector{N, T}`: The operating wavelength(s). For single-wavelength pyrometers (`N=1`), it holds the target wavelength. For spectral band pyrometers (`N=2`), it defines the boundaries `[λₗ, λᵣ]`.
- `ϵ::Base.RefValue{T}`: In-place mutable effective emissivity of the measured target surface.

# Constructors

    Pyrometer(type::Symbol, D::DataType=Float64)

Construct a predefined pyrometer configuration extracted from the `DefaultPyrometersTypes` dictionary.

    Pyrometer(λ::NTuple{N, T}; type::Symbol=:def, ϵ::Number=1.0)

Construct a pyrometer using a tuple of wavelengths `λ`, which is automatically converted into a static vector `SVector`.

    Pyrometer(λ::Union{AbstractVector{T}, T}; type::Symbol=:def, ϵ::Number=1.0)

Universal constructor accepting `λ` either as a single scalar number (monochromatic) or as a vector (band boundaries).

# Examples
```julia
# Create a predefined pyrometer configuration from the package dictionary
p_builtin = Pyrometer(:P)

# Monochromatic pyrometer at 0.65 μm
p_single = Pyrometer(0.65, ϵ=0.85)

# Wide band pyrometer covering 2.4 μm to 8.5 μm
p_band = Pyrometer((2.4, 8.5), type=:mid_ir, ϵ=0.33)
```
"""
struct Pyrometer{N , T} <: AbstractPyrometer{N,T} # this type supports methods for radiative pyrometers
    type::Symbol
    λ::SVector{N,T}
    ϵ::Base.RefValue{T}
            """
        Pyrometer(type::String)

    Pyrometer object Constructor, the type of pyrometer can chosen from the DefaultPyrometersTypes dictionary
    Input:
    type - pyrometer type, must be member of DefaultPyrometersTypes 
    """
    Pyrometer(type::Symbol , D::DataType = Float64) = begin
        if haskey(DefaultPyrometersTypes,type) 
            N = length(DefaultPyrometersTypes[type])
            return new{N , D}(type,
                        DefaultPyrometersTypes[type],
                        Ref(1.0)
                        ) 
        else
                error("Unknown pyrometer type")
        end
    end
    Pyrometer(λ::NTuple{N , T}; type::Symbol=:def ,  ϵ::Number=1.0 ) where {N,T} = begin 
        return new{N , T}(type , SVector{N}(λ)
                    , Ref(ϵ)
                    )
    end
    Pyrometer(λ::Union{AbstractVector{T} , T} ; type::Symbol=:def ,  ϵ::Number=1.0) where T <: Number = begin
        #@assert 0  < ϵ <= 1.0 "Emissivity should be within the (0..1] interval"
        if λ isa AbstractVector
            N = length(λ)
            @assert N == 1 || N == 2 "λ should be a vector of two  Floats or a single Float number"
        else
            N = 1
        end
        new{N , T}(type,SVector{N}(λ),
                                    Ref(ϵ)
            )
    end
end

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


include("MultiwavelengthPyrometryTypes.jl")