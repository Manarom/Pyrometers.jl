
module RadiationPyrometers

    using   LinearAlgebra,
            StaticArrays,
            OrderedCollections, 
            Roots,
            QuadGK

    import  PlanckFunctions as Planck
    export Pyrometer,
        DefaultPyrometersTypes,
        fit_ϵ!,
        fit_ϵ_wavelength!
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
            return new{N , T}(type,SVector{N}(λ)
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
    SingleWavelengthPyrometer

Type alias for a brightness pyrometer operating at **a single, fixed wavelength** (`Pyrometer{1}`). 
Evaluates signals utilizing monochromatic Planck intensity functions (`ibb`).
"""
const SingleWavelengthPyrometer = Pyrometer{1}
"""
    SpectralBandPyrometer

Type alias for a brightness pyrometer operating over **a defined spectral band** (`Pyrometer{2}`). 
Evaluates signals by integrating the Planck distribution function over the defined range (`band_power`).
"""
const SpectralBandPyrometer = Pyrometer{2}
#Constructors
"""
    SingleWavelengthPyrometer(λ::Number; type::Symbol=:def, ϵ::Number=1.0)


"""
SingleWavelengthPyrometer(λ::Number; type::Symbol=:def, ϵ::Number=1.0) = Pyrometer(λ; type=type, ϵ=ϵ)
"""
    SpectralBandPyrometer(λₗ::Number, λᵣ::Number; type::Symbol=:def, ϵ::Number=1.0)


"""
SpectralBandPyrometer(λₗ::Number, λᵣ::Number; type::Symbol=:def, ϵ::Number=1.0) = Pyrometer((λₗ, λᵣ); type=type, ϵ=ϵ)
"""
    SpectralBandPyrometer(λ::NTuple{2}; type::Symbol=:def, ϵ::Number=1.0)


"""
SpectralBandPyrometer(λ::NTuple{2}; type::Symbol=:def, ϵ::Number=1.0) = Pyrometer(λ ; type=type, ϵ=ϵ)
"""
    TwoWavelengthRatioPyrometer

Type alias for a ratio pyrometer measuring the quotient of signals across **two discrete wavelengths**. 
The underlying spectral coordinate type `DT` is expected to be a subtype of `Number`.
"""
const TwoWavelengthRatioPyrometer = RatioPyrometer{2 , T , DT} where {T , DT <: Number}
"""
    TwoBandsRatioPyrometer

Type alias for a ratio pyrometer measuring the quotient of signals across **two distinct spectral bands**. 
The underlying spectral coordinate type `DT` is expected to be a subtype of `Tuple` containing pairs of band edges.
"""
const TwoBandsRatioPyrometer = RatioPyrometer{2 , T , DT} where {T , DT <: Tuple}

# Constructor for TwoWavelengthRatioPyrometer (discrete wavelengths)
function TwoWavelengthRatioPyrometer(λ1::Number, λ2::Number; type::Symbol=:def, ϵ1::Number=1.0, ϵ2::Number=1.0)
    return RatioPyrometer((λ1, λ2); type=type, ϵ1=ϵ1, ϵ2=ϵ2)
end

# Constructor for TwoBandsRatioPyrometer (integrated spectral bands)
function TwoBandsRatioPyrometer(band1::NTuple{2, <:Number}, band2::NTuple{2, <:Number}; 
                               type::Symbol=:def, ϵ1::Number=1.0, ϵ2::Number=1.0)

    return RatioPyrometer((band1, band2); type=type, ϵ1=ϵ1, ϵ2=ϵ2)
end
    """
    wlength(::Pyrometer{N}) where N

Returns the number of wavelengths
"""
wlength(::Pyrometer{N}) where N = N
is_spectral_band(::AbstractPyrometer)  = false
    """
    is_spectral_band(p::Pyrometer)

True if pyrometer `p` is a spectral-band pyrometer (works on a fixed wavelengh region)
"""
is_spectral_band(::SpectralBandPyrometer)  = true
is_spectral_band(::TwoBandsRatioPyrometer) = true
"""
    is_single_wavelength(::Pyrometer{N}) where N

True if Pyrometer is a single wavelength
"""
is_single_wavelength(::SingleWavelengthPyrometer) = true
is_single_wavelength(::TwoWavelengthRatioPyrometer)  = true
is_single_wavelength(::AbstractPyrometer) = false
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
(optional)
`T_starting`  - starting temperature value

"""
function measure(p::AbstractPyrometer , i::D ; T_starting::T=600.0) where {D <: Number, T <: Number}
        ϵ = _get_epsilon_equivalent(p)
        λ = p.λ
        return Roots.find_zero(t -> _Dₜpyro(λ , i , t , ϵ) , T_starting ,  Roots.Halley())  
    end

    Dₜpyro(p::AbstractPyrometer , i , t) = _Dₜpyro(p.λ , i , t , _get_epsilon_equivalent(p))

    _get_epsilon_equivalent(p::Pyrometer)  = p.ϵ[]
    _get_epsilon_equivalent(p::RatioPyrometer)  = e_slope(p)
    (p::AbstractPyrometer)(i; T_starting::Number=1000.0) = measure(p , i , T_starting = T_starting)
    # radiation pyrometer in band 
    _Dₜpyro(λ::SVector{2} , i , t  , ϵ) = _to_halley(Planck.Dₜband_power(t , λₗ = λ[1] , λᵣ = λ[2])  , i , ϵ) 
    # radiation pyrometry for single wavelength
    _Dₜpyro(λ::SVector{1} , i , t  , ϵ) = _to_halley(Planck.Dₜibb(λ[] , t)  , i , ϵ) 
    # spectral ratio pyrometers (single wavelengh)
    _Dₜpyro(λ::NTuple{2 , T} , i , t  , e_slope) where T <: Number = _to_halley(Planck.Dₜspectral_ratio(λ[1] , λ[2] , t , e_slope = 1.0)  , i , e_slope) 
    # spectral ratio band pyrometer 
     _Dₜpyro(λ::NTuple{2 , T} , i , t  , e_slope) where T <: Tuple = _to_halley(Planck.Dₜspectral_band_ratio(λ[1] , λ[2] , t , e_slope = 1.0)  , i , e_slope) 
    """
    _to_halley(tpl , i , ϵ)

Internal function which converts arguments to Roots.jl Halley method from PlanckFunctions Dₜ ... function 
"""
_to_halley(tpl , i , ϵ) = begin 
        (bp , bpd , bpdd) = (tpl[1] , tpl[2] , tpl[3])
        iim = (ϵ *bp - i)
        return ( iim ,  iim / (ϵ * bpd) , bpd/bpdd)
    end

"""
    corrected_temperature(p::Pyrometer, Tmeasured::Number, Tenv::Number, ϵ_env::Number=1.0)

Calculate the actual surface temperature by excluding reflected ambient radiation.

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
- `Tenv::Number`: The temperature of the external radiation source or ambient environment (Kelvins).
- `ϵ_env::Number`: The emissivity of the external radiation source (defaults to 1.0 for a blackbody environment).
"""
function corrected_temperature(p::Pyrometer, Tmeasured::T, Tenv::Number , ϵ_env::Number ) where {T}
        ϵ = _get_epsilon_equivalent(p)
        measured_signal = signal(p , Tmeasured)
        reflected_signal = ϵ_env * (one(T) - ϵ) * signal(p , Tenv) / ϵ
        return p(measured_signal - reflected_signal) 
    end
function corrected_temperature(p::Pyrometer, Tmeasured::T, incident_radiation::Number) where {T}
        ϵ = _get_epsilon_equivalent(p)
        measured_signal = signal(p , Tmeasured)
        reflected_signal =  (one(T) - ϵ) * incident_radiation
        return p(measured_signal - reflected_signal) 
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
    integrate_intensity(p::AbstractPyrometer , λ::AbstractVector , intensity::AbstractVector)

Integral `intensity` signal from discrete data provided for wavelengths `λ`
pyrometer specifies the spectral range 
"""
function integrate_intensity(p::AbstractPyrometer , λ::AbstractVector , intensity::AbstractVector) end
function integrate_intensity(p::SpectralBandPyrometer , λ::AbstractVector , intensity::AbstractVector)
    @assert issorted(λ) "Wavelength vector must be sorted"
    @assert length(λ) == length(intensity) "Vectors must be of the same length"
    (l , f)  = extract_pyrometer_range(p , λ)
    _λ = @view λ[l:f]
    _i = @view intensity[l:f]
    _trapz(_λ , _i)
end

integrate_intensity(p::SingleWavelengthPyrometer ,
                 λ::AbstractVector , intensity::AbstractVector) = _local_interpolate(p.λ[] , λ , intensity)
                 
integrate_intensity(p::TwoWavelengthRatioPyrometer ,
                 λ::AbstractVector , intensity::AbstractVector) = _local_interpolate(p.λ[1] , λ , intensity)/_local_interpolate(p.λ[2] , λ , intensity)

function integrate_intensity(p::TwoBandsRatioPyrometer , λ , intensity) 
    @assert issorted(λ) "Wavelength vector must be sorted"
    @assert length(λ) == length(intensity) "Vectors must be of the same length"
    (i1 , i2)  = ntuple(2) do ii
        (l , f)  = extract_subrange_inds(p.λ[ii]... , λ)
        _λ = @view λ[l:f]
        _i = @view intensity[l:f]
        _trapz(_λ , _i)
    end
    return i1/i2
end
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
            append!(λ,l[2])
            push!(pyr_names,l[1])
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
function full_wavelength_range(p::Vector{Pyrometer})
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
        pyr_vec = Vector{Pyrometer}()
        for (k ,l) in DefaultPyrometersTypes
            push!(pyr_vec,Pyrometer(l[1]  , type = k))
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
    return _integral_emissivity(p.λ[1] , p.λ[2] , λ , ϵ , Tref)
end

function integral_emissivity(p::Union{SpectralBandPyrometer , TwoBandsRatioPyrometer} ,  ϵ_func , Tref::Number)
    return _integral_emissivity(p.λ[1] , p.λ[2] ,  ϵ_func , Tref)
end

integral_emissivity(p::SingleWavelengthPyrometer ,  λ::AbstractVector , ϵ::AbstractVector , _::Number) = (_local_interpolate(p.λ[] , λ , ϵ) , )
integral_emissivity(p::SingleWavelengthPyrometer , ϵ_func , _::Number) = (ϵ_func(p.λ[]) , )
integral_emissivity(p::TwoWavelengthRatioPyrometer ,  λ::AbstractVector , ϵ::AbstractVector , _::Number) = (_local_interpolate(p.λ[1] , λ , ϵ) , _local_interpolate(p.λ[2] , λ , ϵ))
integral_emissivity(p::TwoWavelengthRatioPyrometer , ϵ_func , _::Number)= ϵ_func.(p.λ)

_integral_emissivity(λ1::Number , λ2::Number ,  λ::AbstractVector , ϵ::AbstractVector , Tref) = _integral_emissivity((λ1,) , (λ2,) ,  λ , ϵ , Tref)
function _integral_emissivity(λ1::NTuple{N}, λ2::NTuple{N} ,  λ::AbstractVector , ϵ::AbstractVector , Tref) where N
    ntuple(N) do i
        (f , l) = extract_subrange_inds(λ1[i] , λ2[i] , λ)
        if isnothing(f) || isnothing(l) || f > l
            error(" λ range does not include pyrometer range [$(λ1[i]), $(λ2[i])]")
        end
        _e = @view ϵ[f:l]
        _l = @view λ[f:l]
        return Planck.planck_averaged(_e , _l , Tref)
    end
end

# this version to work with the esmissivity as a function/interpoaltion/polynomial
_integral_emissivity(λ1::Number, λ2::Number  , ϵ_func , Tref) = _integral_emissivity((λ1,), (λ2,)  , ϵ_func , Tref)
function _integral_emissivity(λ1::NTuple{N}, λ2::NTuple{N}  , ϵ_func , Tref) where N
    ntuple(N) do i
        l , r = λ1[i] , λ2[i]
        (emin , emax) = ϵ_func(l) , ϵ_func(r)
        ϵ_baseline = (emin + emax) / 2 # the idea is to exclude the average value 
        irem, _ = quadgk(λ -> (ϵ_func(λ) - ϵ_baseline) * Planck.ibb(λ, Tref), l, r)
        denom = Planck.band_power(Tref , λₗ = l, λᵣ = r )
        return (ϵ_baseline * denom + irem)/denom
    end
end

extract_pyrometer_range(p::SpectralBandPyrometer , λ::AbstractVector) =    extract_subrange_inds(p.λ[1] , p.λ[2] , λ)

extract_subrange_inds(l1 , l2 , λ) = (searchsortedfirst(λ , l1 ) , searchsortedlast( λ , l2))
    """
    fit_ϵ!(p::Vector{Pyrometer},Treal::Float64,Tmeasured::Vector{Float64})

Fits the emissivity of pyrometers to make measured temperature `Tmeasured` fit
fit the real temperature `Treal`

Input:
p - pyrometer objects vector , [Nx0]
Treal - real temperature of the surface, Kelvins
Tmeasured - temperatures measured by the pyrometers, in Kelvins, [Nx0]

"""
function fit_ϵ!(p::Vector{Pyrometer} , Treal::D , Tmeasured::Vector{T}) where {D <: Number ,T <: Number}
        @assert length(p)==length(Tmeasured)  "Vectors must be of the same size"
        N = length(p)
        e_out = Vector{T}(undef , N)
        Threads.@threads for i in 1:N
            @inbounds e_out[i] = fit_ϵ!(p[i] , Tmeasured[i] , Treal)
        end
        return e_out
    end

"""
    fit_ϵ_wavelength!(p::Vector{Pyrometer},Treal::Float64,Tmeasured::Vector{Float64})

The same as [`fit_ϵ!`](@ref) except that it returns the vector of fitted emissivities 
of the same length to the total number of wavelength in all pyrometers in vaector `p`,
e.g. if p[i] is the narrow-band pyrometer 
"""
function fit_ϵ_wavelength!(p::Vector{Pyrometer},Treal::Float64,Tmeasured::Vector{Float64})  
        total_wavelength_number =  sum(wlength,p)
        e_out= Vector{Float64}(undef,total_wavelength_number)
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
function fit_ϵ_wavelength!(p::Pyrometer , Tmeasured::Float64 , Treal::Float64) # this is the same as fit_ϵ! with the exception that 
        # this function returns a vector of values, if pyrometer is single wavelength it returns one -element array
        e_out = similar(p.λ)
        e_out .=fit_ϵ!(p,Tmeasured,Treal)
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
                !isapprox(λ,λp[1],atol=0.05) ? continue : return k
            else
               !(λp[1]<=λ<=λp[2]) ? continue : return k
            end
        end 
        return ""
    end
    Base.show(io::IO, p::Pyrometer{1}) = print(io, "$(p.type) - type: single-wavelength pyrometer:λ = $(p.λ[1]) μm,ϵ = $(p.ϵ[])")
    Base.show(io::IO, p::Pyrometer{2}) = print(io, "$(p.type) - type: spectral-band pyrometer:λ ∈ $(p.λ[1]) ... $(p.λ[2]) μm,ϵ = $(p.ϵ[])")
    Base.show(io::IO, p::RatioPyrometer{2 , T , DT}) where {T , DT <: Number} = print(io, "$(p.type) - type: two wavelength ratio pyrometer:λ₁= $(p.λ[1]) , λ₂ = $(p.λ[2]) μm, ϵ₁ = $(p.ϵ1[]) , ϵ₂ = $(p.ϵ2[]) , e_slope = $(e_slope(p))")
    Base.show(io::IO, p::RatioPyrometer{2 , T , DT}) where {T , DT <: Tuple} = print(io, "$(p.type) - type: two bands ratio pyrometer:λ₁= $(p.λ[1]) , λ₂ = $(p.λ[2]) μm, ϵ₁ = $(p.ϵ1[]) , ϵ₂ = $(p.ϵ2[]) , e_slope = $(e_slope(p))")
    
    
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
# internal function for trapz integration 
    function _trapz(x::AbstractVector, y::AbstractVector)
        N = length(x)
        @assert length(y) == N "Vectors must have the same length"
        N < 2 && return zero(eltype(y))

        s = zero(promote_type(eltype(x), eltype(y), Float64))
        
        @inbounds @simd  for i in 1:(N - 1)
            Δx = x[i+1] - x[i]
            ∑y = y[i+1] + y[i]
            s += Δx * ∑y
        end

        return s * 0.5
    end
end