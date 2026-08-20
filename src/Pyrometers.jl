
module Pyrometers

    using   LinearAlgebra,
            StaticArrays,
            OrderedCollections, 
            Roots,
            QuadGK,
            ADTypes

    import  PlanckFunctions as Planck
    
    export SpectralBandPyrometer, 
        SingleWavelengthPyrometer , 
        TwoBandsRatioPyrometer ,
        TwoWavelengthRatioPyrometer , 
        convert_temperature,
        corrected_temperature,
        integral_emissivity,
        DefaultPyrometersTypes,
        fit_ϵ! , fit_ϵ , 
        Pyrometer , RatioPyrometer , 
        TabularQuantity , AnalyticalSpectralQuantity ,
        IsothermalSpectralQuantity , GenericDifferentiableSpectralQuantity
    """
    Default pyrometers types 

"""
const DefaultPyrometersTypes = OrderedDict(
                    :P => SVector{2}([2.0; 2.6]),
                    :M => SVector{1}([3.4]), 
                    :D => SVector{1}([3.9]),
                    :L => SVector{1}([4.6]),
                    :E => SVector{2}([4.8, 5.2]),
                    :F => SVector{1}([7.9]),
                    :K => SVector{2}([8.0, 9.0]),
                    :B => SVector{2}([9.1,14.0])
    )

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

    abstract type AbstractQuantity{LT , ET} end
"""
    TabularQuantity{LT <: AbstractVector , ET <: AbstractVector} <: AbstractQuantity{LT , ET }

Type wrapper around discrete 

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
    struct TabularQuantity{LT <: AbstractVector , ET <: AbstractVector} <: AbstractQuantity{LT , ET }
        λ::LT
        i::ET
        TabularQuantity(l::LT , e::ET) where {LT <: AbstractVector , ET <: AbstractVector} = begin 
            @assert issorted(l) "First argument must be sorted in ascending order"
            @assert length(l) == length(e) "Two vectors must be of the same length"
            return new{LT , ET}(l , e)
        end
    end
    const IsothermalSpectralQuantity = Planck.IsothermalSpectralQuantity 
    const AnalyticalSpectralQuantity = Planck.AnalyticalSpectralQuantity
    const AbstractSpectralQuantity =  Planck.AbstractSpectralQuantity
    """
        GenericDifferentiableSpectralQuantity(f, backend=AutoForwardDiff())

    An AD-backend agnostic wrapper for a temperature- and wavelength-dependent 
    spectral quantity `f(λ, T)`.

    # Arguments
    * `f`: The core function `f(λ, T)`.
    * `backend`: An `ADTypes.AbstractADType` token specifying the preferred 
    AD framework (e.g., `AutoForwardDiff()`, `AutoZygote()`, `AutoEnzyme()`).
    """
    struct GenericDifferentiableSpectralQuantity{F, B<:AbstractADType} <: AbstractSpectralQuantity
        f::F
        backend::B
        GenericDifferentiableSpectralQuantity(f::F , ad_backend_type::B = AutoForwardDiff()) where {F , B <: AbstractADType} = new{F , B}(f , ad_backend_type)
    end
    """
    Planck.eval_Dₜ(::GenericDifferentiableSpectralQuantity{F, B}, _, _) where {F, B}

Default implementation of GenericDifferentiableSpectralQuantity returns an error
"""
Planck.eval_Dₜ(::GenericDifferentiableSpectralQuantity{F, B} , _ , _) where {F, B} = error("The AD backend $(B) is not loaded in your current session. Please execute `using $(string(B)[5:end])` to activate it.")

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
    measure(p::AbstractPyrometer , i::D ; T_starting::T=600.0) where {D <: Number, T <: Number}

Calculates the "measured" temperature from "measured" signal `i`.
The signal units should be:
$(Planck.units(Planck.ibb)) - for a single wavelength pyrometer , 
$(Planck.units(Planck.band_power)) - for  wide - band pyrometer ,
`dimentionless`  - for sigle wavelength spectral ratio and wide-band spectral ratio; 

# Arguments:
`p` - pyrometer object
`i` - measured signal
`T`  - temperature

"""
function measure(p::AbstractPyrometer , i::D , T::DT=600.0) where {D <: Number, DT <: Number}
        ϵ = _get_epsilon_equivalent(p)
        λ = p.λ
        return Roots.find_zero(t -> _Dₜpyro(λ , i , t , ϵ) , T ,  Roots.Halley())  
    end

    Dₜpyro(p::AbstractPyrometer , i , t) = _Dₜpyro(p.λ , i , t , _get_epsilon_equivalent(p))

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
    get_emissivity(p::RatioPyrometer) = (p.ϵ1[] , p.ϵ[])
    
    (p::AbstractPyrometer)(i; T_starting::Number=1000.0) = measure(p , i , T_starting = T_starting)

    # radiation pyrometer in band 
    _Dₜpyro(λ::SVector{2} , i , t  , ϵ) = _to_halley(Planck.Dₜband_power(t , λₗ = λ[1] , λᵣ = λ[2])  , i , ϵ) 
    # radiation pyrometry for single wavelength
    _Dₜpyro(λ::SVector{1} , i , t  , ϵ) = _to_halley(Planck.Dₜibb(λ[] , t)  , i , ϵ) 
    # spectral ratio pyrometers (single wavelengh)
    _Dₜpyro(λ::NTuple{2 , T} , i , t  , e_slope) where T <: Number = _to_halley(Planck.Dₜspectral_ratio(λ[1] , λ[2] , t )  , i , e_slope) 
    # spectral ratio band pyrometer 
    _Dₜpyro(λ::NTuple{2 , T} , i , t  , e_slope) where T <: Tuple = _to_halley(Planck.Dₜspectral_band_ratio(λ[1] , λ[2] , t )  , i , e_slope) 
    """
    _to_halley(tpl , i , ϵ)

Internal function which converts arguments to Roots.jl Halley method from PlanckFunctions Dₜ ... function 
"""
_to_halley(tpl , i , ϵ) = begin 
        (bp , bpd , bpdd) = (tpl[1] , tpl[2] , tpl[3])
        iim = (ϵ *bp - i)
        return ( iim ,  iim / (ϵ * bpd) , bpd/bpdd)
    end
_to_halley(tpl , i ) = begin 
        (bp , bpd , bpdd) = (tpl[1] , tpl[2] , tpl[3])
        iim = (bp - i)
        return ( iim ,  iim / bpd , bpd/bpdd)
    end    
"""
    subrange_view(λ1::Number , λ2::Number , λ::AbstractVector , i::AbstractVector)

Returns view of two vectors based on `λ[ λ1 <= λ <= λ2]` , `λ` must be sorted in ascending order
"""
function subrange_view(λ1::Number , λ2::Number , e::TabularQuantity)
    (f , l) = extract_subrange_inds(λ1 , λ2 , e.λ)
    if λ1 < first(e.λ) || λ2 > last(e.λ) || f > l
        error("λ range [$(λ1), $(λ2)] goes outside available tabular data bounds [$(first(e.λ)), $(last(e.λ))]")
    end
    _i = @view e.i[f:l]
    _l = @view e.λ[f:l]
    return (_l , _i)
end
subrange_view(λ1::NTuple{2} , λ2::NTuple{2} , e::TabularQuantity) = begin 
    (l1 , e1) = subrange_view(λ1[1] , λ1[2] , e::TabularQuantity) 
    (l2 , e2) = subrange_view(λ2[1] , λ2[2] , e::TabularQuantity) 
    ( (l1 , l2)
        , 
       (e1 , e2)
    )
end
extract_pyrometer_inds(p::SpectralBandPyrometer , λ::AbstractVector ) =    extract_subrange_inds(p.λ[1] , p.λ[2] , λ)
extract_subrange_inds(l1 , l2 , λ) = (searchsortedfirst(λ , l1 ) , searchsortedlast( λ , l2))
"""
Internal object, which wrapps the measured signal
"""
abstract type AbstractQuantityContext{L , E , F} end
"""
    Type wrapper for [`TabularQuantity`](@ref) which wrapps the measured spectrum
"""
struct TabularQuantityContext{ L , E , F} <: AbstractQuantityContext{L , E , F}
    _l::L
    _e::E
    i_measured::F
# main constructor l1 - number for 
TabularQuantityContext(l1 , l2 , e::TabularQuantity , imeasured::F) where F = begin 
        (_l , _e) = subrange_view(l1  , l2 , e)
        new{typeof(_l) , typeof(_e) , F}(_l , _e, imeasured)
    end
    """
    TabularQuantityContext(l::NTuple{2 , D} , e::TabularQuantity , imeasured) where D <: Tuple

version of constructor with input wavelength for two-band ratio pyrometer 
"""
TabularQuantityContext(l::NTuple{2 , D} , e::TabularQuantity , imeasured) where D <: Tuple = begin 
            TabularQuantityContext(l[1] , l[2] , e , imeasured)
    end
end
"""
    TabularQuantityContext(l::StaticArray{Tuple{2}, D, 1},e::TabularQuantity , imeasured) where D <: Number

Additional constructor for spectral-band pyrometer, which stores wavelength as `SVector`
"""
TabularQuantityContext(l::StaticArray{Tuple{2}, D, 1},e::TabularQuantity , imeasured) where D <: Number  = TabularQuantityContext(l[1] , l[2] , e , imeasured) 


(ctx::TabularQuantityContext)(t) = _to_halley(Planck.Dₜplanck_weighted(ctx._e , ctx._l , t) , ctx.i_measured)
# this version of wrapper is for the TwoBandsRatioPyrometer
function (ctx::TabularQuantityContext{L , E})(t) where {L <: Tuple , E <: Tuple}
    tpl = Planck.Dₜplanck_weighted_ratio(ctx._e[1] , ctx._l[1] ,ctx._e[2] , ctx._l[2] , t )
    return _to_halley(tpl, ctx.i_measured )
end

"""
    measure(p::AbstractPyrometer , imeasured::Number , ϵ::AbstractSpectralQuantity; T_starting::Number = 600.0)

    General function to measure the temperature from the signal `imeasured` taking into account 
the emissivity provided as [`TabularQuantity`](@ref), [``]
    
# Arguments
- `p`: AbstractPyrometer object 


# Examples
```julia
# Spectral band pyrometer for 2.0 - 4.5 μm spectral band 
p_band = SpectralBandPyrometer((2.0 , 4.5) , (7.0 , 9.1) , ϵ1 = 0.4 , ϵ2 = 0.93)
```
"""
measure(p::AbstractPyrometer , imeasured::Number , ϵ::AbstractSpectralQuantity , T::Number = 600.0) = error("not implemented yet")
 
function measure(p::Union{SpectralBandPyrometer , TwoBandsRatioPyrometer} , imeasured::Number , ϵ::TabularQuantity, T::Number = 600.0)
    ctx = TabularQuantityContext(p.λ , ϵ , imeasured)
    return Roots.find_zero(ctx , T ,  Roots.Halley())
end 
function measure(p::Union{SingleWavelengthPyrometer , TwoWavelengthRatioPyrometer} ,
                imeasured::Number , ϵ::AbstractSpectralQuantity ,  
                t::Number = 600.0)

    e_previous = get_emissivity(p)
    e_new = _local_interpolate(p.λ[] , ϵ.λ , ϵ.i)
    set_emissivity!(p , e_new)
    val = p(imeasured)
    set_emissivity!(p , e_previous)
    return val
end
_get_single_wavelength_value(l, _ , e::TabularQuantity) = _local_interpolate(l , e.λ , e.i)
_get_single_wavelength_value(l::Number, _ , e::TabularQuantity) = _local_interpolate(l , e.λ , e.i)
function measure(p::Union{SpectralBandPyrometer , TwoBandsRatioPyrometer} , 
                    imeasured::Number , ϵ::AbstractSpectralQuantity; 
                    T_starting::Number = 600.0)

    ctx = GenericSpectralQuantityContext(Tuple(p.λ) , ϵ , imeasured)
    return Roots.find_zero(ctx , T_starting ,  Roots.Halley())
end 

struct GenericSpectralQuantityContext{L , E , F} <: AbstractQuantityContext{L , E , F}
    λ::L
    e::E
    i_measured::F
    GenericSpectralQuantityContext(λ::L , quantity::E , i_measured::F) where {L , E , F} = new{L , E , F}(λ , quantity , i_measured)
end

(ctx::AbstractQuantityContext{L})(t) where {L <: NTuple{2 , D}} where D <: Number = _to_halley(Planck.Dₜplanck_weighted(ctx.e ,ctx.λ[1] ,ctx.λ[2], t) , ctx.i_measured)
# this version of wrapper is for the TwoBandsRatioPyrometer
function (ctx::AbstractQuantityContext{L})(t) where {L <: NTuple{2 , D} } where D <: NTuple 
    tpl = Planck.Dₜplanck_weighted_ratio(ctx.e , ctx.λ[1] , ctx.λ[2] , t )
    return _to_halley(tpl, ctx.i_measured )
end










"""
    stray_radiation_corrected_temperature(p::Pyrometer, Tmeasured::T; Tenv::Number , ϵ_env::Number ) where {T}

Calculate the actual surface temperature by excluding stray radiation reflected to the pyrometers's FOV.

If measurements are performed in an environment with a background temperature `Tenv`
and a source emissivity `ϵ_env`, the pyrometer will overestimate the true surface
temperature, reporting `Tmeasured` instead. 

**supposes that the pyrometer emissivity is settled to the emissivity of the surface**

This function:
1. Converts both `Tmeasured` and `Tenv` to corresponding radiation signals 
   (integral band power or spectral intensity, depending on the pyrometer type).
2. Subtracts the reflected portion of the ambient radiation from the total measured signal.
3. Recalculates and returns the true, reflection-free surface temperature.

# Arguments
- `p::Pyrometer`: The pyrometer object containing surface emissivity and wavelength parameters.
- `Tmeasured::Number`: The raw, uncorrected temperature reading from the pyrometer (Kelvins).
(kwargs)
- `Tenv::Number`: The temperature of the external radiation source or ambient environment (Kelvins).
- `ϵ_env::Number`: The emissivity of the external radiation source (defaults to 1.0 for a blackbody environment).
"""
function stray_radiation_corrected_temperature(p::Pyrometer, Tmeasured::T; Tenv::Number , ϵ_env::Number ) where {T}
        return stray_radiation_corrected_temperature(p , Tmeasured , ϵ_env * signal(p , Tenv) / _get_epsilon_equivalent(p))
    end
"""
    stray_radiation_corrected_temperature(p::Pyrometer, Tmeasured::T, incident_radiation::Number) where {T}

    The  as  [`stray_radiation_corrected_temperature(p::Pyrometer, Tmeasured::T; Tenv::Number , ϵ_env::Number )`](@ref)
but external radiation is provided as a single value (it should be integrated within the pyrometer working range) 
"""
function stray_radiation_corrected_temperature(p::Pyrometer, Tmeasured::T, incident_radiation::Number) where {T}
        ϵ = _get_epsilon_equivalent(p)
        measured_signal = signal(p , Tmeasured)
        reflected_signal =  (one(T) - ϵ) * incident_radiation
        return p(measured_signal - reflected_signal) 
    end 

function stray_radiation_corrected_temperature(p::Pyrometer, Tmeasured::T, incident_radiation_function) where {T}
        return stray_radiation_corrected_temperature(p::Pyrometer, Tmeasured::T, integrate_intensity(p ,intensity_function = incident_radiation_function))
    end     
    """
    signal(p::AbstractPyrometer , Tmeasured::Number)

Returns the signal value which will give the temperature `Tmeasured`
"""
function signal(p::AbstractPyrometer , Tmeasured::Number) end
signal(p::SingleWavelengthPyrometer , Tmeasured::Number)  = p.ϵ[] * Planck.ibb(p.λ[] , Tmeasured)
signal(p::SpectralBandPyrometer , Tmeasured::Number) = p.ϵ[] * Planck.band_power(Tmeasured , λₗ=p.λ[1] , λᵣ=p.λ[2])
signal(p::TwoWavelengthRatioPyrometer , Tmeasured::Number)  = Planck.spectral_ratio(p.λ[1] , p.λ[2] , Tmeasured , e_slope = e_slope(p))
signal(p::TwoBandsRatioPyrometer , Tmeasured::Number) = Planck.spectral_band_ratio(p.λ[1] , p.λ[2] , Tmeasured , e_slope = e_slope(p))

"""
    integrate_intensity(p::AbstractPyrometer ; λ , intensity)

Converts the input spectral `intensity` given as a discrete set of values for wavelengths `λ`
to the pyrometer signal.

Returns the quantity, which is equal to the type of pyrometer signal , e.g. if pyrometer is 

- `SingleWavelengthPyrometer` returns the spectral intensity at the wavelength of pyrometer

- `SpectralBandPyrometer`  - total intensity within the spectral range the pyrometer

- `TwoWavelengthRatioPyrometer`  - ratio of intensities at two wavelength 

- `TwoBandsRatioPyrometer`  - ratio of total intensities for two spectral ranges of the pyrometer


`λ` , `intensity`  (both at the same time) if the intensity is provided as a discrete set of points 

"""
function integrate_intensity(p::AbstractPyrometer ; λ::AbstractVector , intensity::AbstractVector) error("This state must not be reached") end
"""
    integrate_intensity(p::AbstractPyrometer ,  intensity_function)

If the intensity is provided as a function of wl
"""
function integrate_intensity(p::AbstractPyrometer ,  intensity_function) error("This state must not be reached") end


integrate_intensity(p::SpectralBandPyrometer ; λ::AbstractVector , intensity::AbstractVector) = _integrate_within_wavelength(p.λ[1] , p.λ[2] , λ , intensity)
integrate_intensity(p::SpectralBandPyrometer , intensity_function) = _integrate_within_wavelength(p.λ[1] , p.λ[2] , intensity_function)

function integrate_intensity(p::TwoBandsRatioPyrometer , intensity_function) 
    (i1 , i2)  = ntuple(2) do ii
        _integrate_within_wavelength(p.λ[ii]... , intensity_function)
    end
    return i1/i2
end
# versions for single wavelength pyrometers 
integrate_intensity(p::SingleWavelengthPyrometer ;
                 λ::AbstractVector , intensity::AbstractVector) = _local_interpolate(p.λ[] , λ , intensity)
integrate_intensity(p::SingleWavelengthPyrometer, intensity_function) = intensity_function(p.λ[])
integrate_intensity(p::TwoWavelengthRatioPyrometer, intensity_function) = intensity_function(p.λ[1])/intensity_function(p.λ[2])
integrate_intensity(p::TwoWavelengthRatioPyrometer ;
                 λ::AbstractVector , intensity::AbstractVector) = _local_interpolate(p.λ[1] , λ , intensity)/_local_interpolate(p.λ[2] , λ , intensity)


"""
    _integrate_within_wavelength(λ1 , λ2 , λ , i)

Simpson's method integration of discrete signal `i` given at wavelengths `λ`
"""
function _integrate_within_wavelength(λ1 , λ2 , λ , intensity)
    @assert issorted(λ) "Wavelength vector must be sorted"
    @assert length(λ) == length(intensity) "Vectors must be of the same length"
    return _simpson(subrange_view(λ1 , λ2 , λ , intensity)...)
end
"""
    _integrate_within_wavelength(λ1 , λ2 , intensity_function)

Using quadgk if the signal is provided as a callable
"""
_integrate_within_wavelength(λ1 , λ2 , intensity_function) = first(quadgk(intensity_function , λ1 , λ2 ))

    """
    fit_ϵ(p::AbstractPyrometer , Tmeasured::Number , Treal::Number)

Finds the emissivity or e_slope 
# Arguments :

`p` - pyrometer object
`Treal` - real temperature of the surface, Kelvins
`Tmeasured` - temperature measured by the pyrometer, Kelvins
"""
function fit_ϵ(p::AbstractPyrometer , Tmeasured::Number , Treal::Number)
         return _get_epsilon_equivalent(p) * signal(p , Tmeasured)/signal(p , Treal)
    end
fit_ϵ!(p::AbstractPyrometer , Tmeasured::Number , Treal::Number) = set_emissivity!(p , fit_ϵ(p , Tmeasured , Treal))
    """
    convert_temperature(p::AbstractPyrometer , Tmeasured  , ϵ_new)

Converts temperature `Tmeasured` measured using pyrometer `p` with it specified emissivity 
to a new temperature measured with `ϵ_new` , the type of `ϵ_new` depends on the type of pyrometer 
if `ϵ_new` is a `Number` than if p is `RatioPyrometer` it assumes `e_new` is `e_slope`, if 
`e_new` is `NTuple{2 , Number}` it modifies both emissivities at two wavelength
"""
function convert_temperature(p::AbstractPyrometer , Tmeasured  , ϵ_new)
    i = signal(p , Tmeasured)
    _e = _get_epsilon_equivalent(p)
    set_emissivity!(p , ϵ_new)
    Tnew = measure(p , i)
    set_emissivity!(p , _e) # returning previous emissivity
    return Tnew
end
    """
    Base.isless(p1::Pyrometer,p2::Pyrometer)

Vector of Pyrometer objects can be sorted using isless
"""
function Base.isless(p1::Pyrometer,p2::Pyrometer) # is used to sort the vector of pyrometers
        return all(p1.λ .< p2.λ)
    end
    """
    wavelength_number()

Returns the length of wavelengths vector covered by the [`DefaultPyrometersTypes`](@ref) (all default pyrometers wavelengh region)
"""
function wavelengths_number()
        return mapreduce(x->length(x),+,DefaultPyrometersTypes)
    end
"""
    wavelengths_number(p::Vector{Pyrometer})

Returns the total number of wavelength for the vector of pyrometers
"""
function wavelengths_number(p::Base.AbstractVecOrTuple{D}) where D <: AbstractPyrometer
    return sum(wlength , p)
end
    """
    full_wavelength_range()

Creates the wavelengths vector covered by default pyrometers see [`DefaultPyrometersTypes`](@ref) 
"""
function full_wavelength_range()
        #sz = mapreduce(x->length(x),+,DefaultPyrometersTypes)
        λ = Vector{Float64}()
        pyr_names = Vector{String}()
        for l in DefaultPyrometersTypes
            append!(λ , l[2])
            push!(pyr_names , l[1])
            if length(l[2])>1
                push!(pyr_names,l[1])
            end
        end
        inds = sortperm(λ)
        pyr_names.=pyr_names[inds]
        λ.=λ[inds]
        return λ,pyr_names
    end
    """
    full_wavelength_range()

Creates the wavelengths vector covered by all pyrometers in vector `p`
"""
function full_wavelength_range(p::Vector{T}) where T <: AbstractPyrometer
        #sz = mapreduce(x->length(x),+,DefaultPyrometersTypes)
        λ = Vector{Float64}(undef,wavelengths_number(p))
        counter = 0
        for pj in p
            if is_spectral_band(pj)
                counter += 2
                λ[counter-1] = pj.λ[1]
            else
                counter+=1
            end    
            λ[counter] = pj.λ[end]
        end
        return λ
    end   
    """
    produce_pyrometers()

Creates the vector of all default pyrometers 
"""
function produce_pyrometers()
        pyr_vec = Vector{AbstractPyrometer}()
        for k in keys(DefaultPyrometersTypes)
            push!(pyr_vec, Pyrometer(k))
        end
        sort!(pyr_vec) # sorting according to the wavelength increase
        return pyr_vec
    end
    """
    set_emissivity!(p::Pyrometer,em_value)

Setter for spectral emissivity
"""
set_emissivity!(p::Pyrometer , em_value::Number) = (p.ϵ[] = em_value)
set_emissivity!(p::Pyrometer , em_value::NTuple{1}) = (p.ϵ[] = em_value[1])
"""
    set_emissivity!(p::RatioPyrometer , em_value::Number)

For spectral ratio pyrometer  `em_value` is the e-slope (emissivities ratio)
"""
set_emissivity!(p::RatioPyrometer , em_value::Number) = begin 
    p.ϵ1[]  = em_value * p.ϵ2[] 
end
"""
    set_emissivity!(p::RatioPyrometer{N , T} , em_value::NTuple{2 , D}) where {D <: Number , N , T}

If emissivity is provided as a tuple it is considered as two emissivities at two spectral ranges of 
ratio pyrometer
"""
set_emissivity!(p::RatioPyrometer{N , T} , em_value::NTuple{2 , D}) where {D <: Number , N , T} = begin 
    p.ϵ1[]  = T(em_value[1])
    p.ϵ2[]  = T(em_value[2]) 
end
"""
    set_emissivity!(p::AbstractPyrometer; λ::AbstractVector , ϵ::AbstractVector , Tref::Number=1273.15)

Setting the emissivity from experimental data of `ϵ` measured at wavelength `λ`

# Arguments
`λ` - wavelength , 
`ϵ` - spectral emissivity 
"""
function set_emissivity!(p::AbstractPyrometer; λ::AbstractVector , ϵ::AbstractVector , Tref::Number=1273.15)
    @assert length(λ) == length(ϵ) " λ::AbstractVector , ϵ::AbstractVector must be of the same length "
    e_int = integral_emissivity(p ,  λ , ϵ , Tref)
    set_emissivity!(p , e_int)
end

"""
    integral_emissivity(p::SpectralBandPyrometer ,  λ::AbstractVector , ϵ::AbstractVector , Tref::Number)

Evaluate the integral emissivity within the pyrometer's working range using external spectral emissivity
provided as discrete values `ϵ` for wavelengths `λ`. 
    
**`λ` vector must be sorted in ascending order**

# Note
This discrete version is highly efficient when the number of wavelength points within the pyrometer's 
working range is small (typically less than 100). For dense spectral datasets (more than 100 points), 
it is faster to interpolate `ϵ` and evaluate the integral using adaptive numerical integration. 
See: [`integral_emissivity(p::Union{SpectralBandPyrometer, TwoBandsRatioPyrometer}, ϵ_func, Tref::Number)`](@ref).

# Returns
`ϵ_int` computed as:
`ϵ_int = ∫ϵ(λ)ibb(λ , Tref)dλ / ∫ibb(λ , Tref)dλ`
"""
function integral_emissivity(p::Union{SpectralBandPyrometer , TwoBandsRatioPyrometer} ,  λ::AbstractVector , ϵ::AbstractVector , Tref::Number)
    @assert issorted(λ) "Wavelengths vector must be sorted" 
    return _pyrometer_band_averaged(p.λ[1] , p.λ[2] , λ , ϵ , Tref)
end

function integral_emissivity(p::Union{SpectralBandPyrometer , TwoBandsRatioPyrometer} ,  ϵ_func , Tref::Number)
    return _pyrometer_band_averaged(p.λ[1] , p.λ[2] ,  ϵ_func , Tref)
end
"""
    integral_temperature_dependent_emissivity(p::Union{SpectralBandPyrometer , TwoBandsRatioPyrometer} , 
     ϵ_func , T::Number)

Evaluate the integral emissivity within the pyrometer's working range using external spectral emissivity
provided as a function of wavelength and temperature `ϵ(λ , T)`. 
    
# Returns
`ϵ_int` computed as:
`ϵ_int(T) = ∫ϵ(λ , T)ibb(λ , T)dλ / ∫ibb(λ , T)dλ`
"""
function integral_temperature_dependent_emissivity(p::Union{SpectralBandPyrometer , TwoBandsRatioPyrometer} , 
     ϵ_func , T::Number)
    return _pyrometer_band_temperature_dependent_averaged(p.λ[1] , p.λ[2] ,  ϵ_func , T)
end


integral_emissivity(p::SingleWavelengthPyrometer ,  λ::AbstractVector , ϵ::AbstractVector , _::Number) = (_local_interpolate(p.λ[] , λ , ϵ) , )
integral_emissivity(p::SingleWavelengthPyrometer , ϵ_func , _::Number) = (ϵ_func(p.λ[]) , )
integral_emissivity(p::TwoWavelengthRatioPyrometer ,  λ::AbstractVector , ϵ::AbstractVector , _::Number) = (_local_interpolate(p.λ[1] , λ , ϵ) , _local_interpolate(p.λ[2] , λ , ϵ))
integral_emissivity(p::TwoWavelengthRatioPyrometer , ϵ_func , _::Number)= ϵ_func.(p.λ)


_pyrometer_band_averaged(λ1::Number , λ2::Number ,  λ::AbstractVector , ϵ::AbstractVector , Tref) = _pyrometer_band_averaged((λ1,) , (λ2,) ,  λ , ϵ , Tref)
_pyrometer_band_temperature_dependent_averaged(λ1::Number , λ2::Number ,  ϵ , Tref) = _pyrometer_band_temperature_dependent_averaged((λ1,) , (λ2,) ,   ϵ , Tref)

function _pyrometer_band_averaged(λ1::NTuple{N}, λ2::NTuple{N} ,  λ::AbstractVector , ϵ::AbstractVector , Tref) where N
    ntuple(N) do i
        ( _l , _e ) = subrange_view(λ1[i] , λ2[i] , λ , ϵ)
        return Planck.planck_averaged(_e , _l , Tref)
    end
end

# this version to work with the esmissivity as a function/interpoaltion/polynomial
_pyrometer_band_averaged(λ1::Number, λ2::Number  , ϵ_func , Tref) = _pyrometer_band_averaged((λ1,), (λ2,)  , ϵ_func , Tref)
function _pyrometer_band_averaged(λ1::NTuple{N}, λ2::NTuple{N}  , ϵ_func , Tref) where N
    ntuple(N) do i
        l , r = λ1[i] , λ2[i]
        (emin , emax) = ϵ_func(l) , ϵ_func(r)
        ϵ_baseline = (emin + emax) / 2 # the idea is to exclude the average value 
        irem, _ = quadgk(λ -> (ϵ_func(λ) - ϵ_baseline) * Planck.ibb(λ, Tref), l, r)
        denom = Planck.band_power(Tref , λₗ = l, λᵣ = r )
        return (ϵ_baseline * denom + irem)/denom
    end
end
"""
    _pyrometer_band_tdependent_averaged(λ1::NTuple{N}, λ2::NTuple{N}  , ϵ_func , T) where N

Averages if the function is both wl and temperature dependent
"""
function _pyrometer_band_temperature_dependent_averaged(λ1::NTuple{N}, λ2::NTuple{N}  , func , T) where N
    ntuple(N) do i
        l , r = λ1[i] , λ2[i]
        (emin , emax) = func(l , T) , func(r , T)
        ϵ_baseline = (emin + emax) / 2 # the idea is to exclude the average value 
        irem, _ = quadgk(λ -> (func(λ , T) - ϵ_baseline) * Planck.ibb(λ, T), l, r)
        denom = Planck.band_power(T , λₗ = l, λᵣ = r )
        return (ϵ_baseline * denom + irem)/denom
    end
end


    """
    fit_ϵ!(p::Vector{P} , Treal::D , Tmeasured::Vector{T}) where {D <: Number ,T <: Number , P <: AbstractPyrometer}


Fits the emissivity of pyrometers to make measured temperature `Tmeasured` fit
fit the real temperature `Treal`

Input:
p - pyrometer objects vector , [Nx0]
Treal - real temperature of the surface, Kelvins
Tmeasured - temperatures measured by the pyrometers, in Kelvins, [Nx0]

"""
function fit_ϵ!(p::Vector{P} , Treal::D , Tmeasured::Vector{T}) where {D <: Number ,T <: Number , P <: AbstractPyrometer}
        @assert length(p)==length(Tmeasured)  "Vectors must be of the same size"
        N = length(p)
        e_out = Vector{T}(undef , N)
        Threads.@threads for i in 1:N
            @inbounds e_out[i] = fit_ϵ!(p[i] , Tmeasured[i] , Treal)
        end
        return e_out
    end

"""
    fit_ϵ_wavelength!(p::Vector{T},Treal::D,Tmeasured::Vector{D})   where {T<: AbstractPyrometer , D <: Number}

The same as [`fit_ϵ!`](@ref) except that it returns the vector of fitted emissivities 
of the same length to the total number of wavelength in all pyrometers in vaector `p`,
e.g. if p[i] is the narrow-band pyrometer 
"""
function fit_ϵ_wavelength!(p::Vector{T},Treal::D,Tmeasured::Vector{D})   where {T<: AbstractPyrometer , D <: Number}
        total_wavelength_number =  sum(wlength,p)
        e_out= Vector{D}(undef,total_wavelength_number)
        counter = 0
        for (i,e) in enumerate(fit_ϵ!(p,Treal,Tmeasured))
            if is_spectral_band(p[i]) 
                counter +=2 
                e_out[counter-1] = e
            else
                counter +=1 
            end
            e_out[counter] = e
        end
        return e_out
    end

    """
    fit_ϵ_wavelength!(p::Pyrometer,Tmeasured::Float64,Treal::Float64)

Fits emissivity and returns it as a vector of the same size as pyrometer's wavelength region
Some pyrometers has 2-wavelength, other work on a single wavelength, for two-wavelength pyrometers
e_out will be a two-element vector
Input:
    p - pyrometer object
    Tmeasured - temperature measured by the pyrometer, Kelvins
    Treal - real temeprature of the surface, Kelvins 
"""
function fit_ϵ_wavelength!(p::AbstractPyrometer , Tmeasured::Float64 , Treal::Float64) # this is the same as fit_ϵ! with the exception that 
        # this function returns a vector of values, if pyrometer is single wavelength it returns one -element array
        e_out = similar(p.λ)
        e_out .= fit_ϵ!(p , Tmeasured , Treal)
        return e_out
    end
    """
    switch_the_type(λ::Float64)

Returns the type (which can be used as a key of DefaultPyrometersTypes dict) depending on the wavelengh
Input:
    λ - wavelength in μm
"""
function switch_the_type(λ::Float64)
        for (k,λp) in DefaultPyrometersTypes
            if length(λp)==1
                !isapprox(λ  , λp[1] , atol=0.05) ? continue : return k
            else
               !(λp[1] <= λ <= λp[2]) ? continue : return k
            end
        end 
        return ""
    end

    Base.show(io::IO, p::SingleWavelengthPyrometer) = print(io, "$(p.type) - type: single-wavelength pyrometer:λ = $(p.λ[1]) μm,ϵ = $(p.ϵ[])")
    Base.show(io::IO, p::SpectralBandPyrometer) = print(io, "$(p.type) - type: spectral-band pyrometer:λ ∈ $(p.λ[1]) ... $(p.λ[2]) μm,ϵ = $(p.ϵ[])")
    Base.show(io::IO, p::TwoWavelengthRatioPyrometer) = print(io, "$(p.type) - type: two wavelength ratio pyrometer:λ₁= $(p.λ[1]) , λ₂ = $(p.λ[2]) μm, ϵ₁ = $(p.ϵ1[]) , ϵ₂ = $(p.ϵ2[]) , e_slope = $(e_slope(p))")
    Base.show(io::IO, p::TwoBandsRatioPyrometer)  = print(io, "$(p.type) - type: two bands ratio pyrometer:λ₁= $(p.λ[1]) , λ₂ = $(p.λ[2]) μm, ϵ₁ = $(p.ϵ1[]) , ϵ₂ = $(p.ϵ2[]) , e_slope = $(e_slope(p))")
    
    _local_interpolate(l::NTuple{N} ,  λ_nodes::AbstractVector, ϵ_nodes::AbstractVector) where N = ntuple(N) do i
         _local_interpolate(l[i] , λ_nodes , ϵ_nodes)
    end
     """
    _local_interpolate(λ_target::Number, λ_nodes::AbstractVector, ϵ_nodes::AbstractVector)

single point interpolation 
"""
@inline function _local_interpolate(λ_target::Number, λ_nodes::AbstractVector, ϵ_nodes::AbstractVector)
       
        idx = searchsortedfirst(λ_nodes, λ_target)
        
        if idx == 1
            return ϵ_nodes[begin]
        elseif idx > length(λ_nodes)
            return ϵ_nodes[end]
        end
        
        λ_start, λ_end = λ_nodes[idx-1], λ_nodes[idx]
        ϵ_start, ϵ_end = ϵ_nodes[idx-1], ϵ_nodes[idx]
        
        t = (λ_target - λ_start) / (λ_end - λ_start)
        return ϵ_start + t * (ϵ_end - ϵ_start)
    end
    include("custom_integration_funcs.jl")
end