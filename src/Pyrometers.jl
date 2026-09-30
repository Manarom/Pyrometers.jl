
module Pyrometers

    using   LinearAlgebra,
            StaticArrays,
            OrderedCollections, 
            Roots,
            QuadGK,
            ADTypes,
            ScaledPolynomials

    import  PlanckFunctions as Planck
   
    export SpectralBandPyrometer, 
        SingleWavelengthPyrometer , 
        TwoBandsRatioPyrometer ,
        TwoWavelengthRatioPyrometer , 
        MultiWavelengthPyrometer , 
        convert_temperature,
        integral_emissivity,
        set_emissivity! , 
        DefaultPyrometersTypes,
        fit_ϵ! , fit_ϵ , 
        Pyrometer , RatioPyrometer , 
        TabularQuantity , AnalyticalSpectralQuantity ,
        IsothermalSpectralQuantity , GenericDifferentiableSpectralQuantity , 
        stray_radiation_corrected_temperature , PlanckEmitter , 
        external_source_corrected_temperature , 
        EnclosureGeometry , 
        ViewFactorGeometry , 
        fix_temperature , 
        SpectralReflectivity,
        emissivity , 
        emissivity_poly
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
                    :B => SVector{2}([9.1 , 14.0])
    )
    abstract type AbstractDiscreteQuantity{LT , ET} end
    const IsothermalSpectralQuantity = Planck.IsothermalSpectralQuantity 
    const AnalyticalSpectralQuantity = Planck.AnalyticalSpectralQuantity
    const AbstractSpectralQuantity =  Planck.AbstractSpectralQuantity

    include("spectral_quantities.jl")
    include("pyrometers_types.jl")
    include("measure_func.jl")
    include("signal_func.jl")
    include("convert_temperature_func.jl")
    include("with_pyrometers_integrator.jl")
    include("stray_radiation_correction_func.jl")   

#basic types functions 


  


"""
    external_source_corrected_temperature(p::AbstractPyrometer, 
                                Tmeasured::Number, 
                                ϵ_surf::Union{IsothermalSpectralQuantity ,  Number},
                                Tsource::Number , 
                                ϵ_src::Union{IsothermalSpectralQuantity ,  Number}  , 
                                ::EnclosureGeometry)

Computes the true surface temperature of an object by isolating and removing parasitic reflected 
radiation originating from an external heated source (e.g., furnace walls, refractory lining) 
within the pyrometer's spectral operating range.
```
"""
function external_source_corrected_temperature(p::AbstractPyrometer, 
                                Tmeasured::Number, 
                                ϵ_surf::Union{IsothermalSpectralQuantity ,  Number},
                                Tsource::Number , 
                                ϵ_src::Union{IsothermalSpectralQuantity ,  Number}  , 
                                ::EnclosureGeometry)

        measured_signal = signal(p, Tmeasured , ϵ_surf) 
        r_eff = ϵ_src * SpectralReflectivity(ϵ_surf)
        reflected_signal = signal(p, Tsource, r_eff)
        pure_signal = _extract_signals(p , measured_signal , reflected_signal)
    
        return p(pure_signal, ϵ_surf)
    end

"""
    external_source_corrected_temperature(p::Pyrometer, Tmeasured::Number, 
                    Tsource::Number ,
                    ϵ_src::Union{Number , AbstractSpectralQuantity}  , 
                    geometry::AbstractRadiationGeometry=EnclosureGeometry())

A convenience forwarding method that extracts the baseline emissivity automatically from the pyrometer (`get_emissivity(p)`) and applies the universal multi-reflection stray radiation correction.

# Arguments
- `p::Pyrometer`: The target single-channel pyrometer instance.
- `Tmeasured::Number`: The raw, uncorrected temperature reading (in Kelvin).
- `Tsource::Number`: The temperature of the parasitic external background source (in Kelvin).
- `ϵ_src::Union{Number, AbstractSpectralQuantity}`: The background source emissivity (constant or functional).
- `geometry::AbstractRadiationGeometry`: Geometric parameter configuration layout (defaults to `EnclosureGeometry()`).
"""
external_source_corrected_temperature(p::Pyrometer, Tmeasured::Number, 
                    Tsource::Number ,
                    ϵ_src::Union{Number , AbstractSpectralQuantity}  = 1.0, 
                    geometry::AbstractRadiationGeometry=EnclosureGeometry()) = external_source_corrected_temperature(p , 
                                                                                            Tmeasured , get_emissivity(p) , 
                                                                                            Tsource , ϵ_src , geometry)



"""
    integral_emissivity(p::AbstractPyrometer , ϵ::AbstractContinuousOrDiscreteQuantity , Tref::Number)

Evaluate the integral emissivity within the pyrometer's working range using external spectral emissivity
provided as discrete or continuous `ϵ` 
    
**`λ` vector must be sorted in ascending order**
# Arguments
`p` - pyrometer object 
`ϵ` - emissivity data
`Tref` - reference temperature
# Returns tuple of single or two numbers depending on pyrometer's type 
E.g.  for [`SpectralBandPyrometer`](@ref) `ϵ_int` is computed as:
`ϵ_int = ∫ϵ(λ)ibb(λ , Tref)dλ / ∫ibb(λ , Tref)dλ`
"""
integral_emissivity(p::AbstractPyrometer , ϵ::AbstractContinuousOrDiscreteQuantity , Tref::Number) = error(" Not specified for p $(typeof(p))")
integral_emissivity(p::Union{SpectralBandPyrometer , TwoBandsRatioPyrometer} ,  
     ϵ::AbstractContinuousOrDiscreteQuantity , Tref::Number) = _pyrometer_band_averaged(p.λ[1] , p.λ[2] , ϵ , Tref)

integral_emissivity(p::Union{SingleWavelengthPyrometer, TwoWavelengthRatioPyrometer} ,   ϵ::AbstractContinuousOrDiscreteQuantity , T_ref::Number) = get_single_wavelength_value(Tuple(p.λ) , T_ref , ϵ)

_pyrometer_band_averaged(λ1::Number , λ2::Number ,   ϵ::AbstractContinuousOrDiscreteQuantity , Tref) = _pyrometer_band_averaged((λ1,) , (λ2,) ,   ϵ , Tref)
function _pyrometer_band_averaged(λ1::NTuple{N}, λ2::NTuple{N} ,   ϵ::AbstractDiscreteQuantity , Tref) where N
    ntuple(N) do i
        ( _l , _e ) = subrange_view(λ1[i] , λ2[i] ,  ϵ)
        return Planck.planck_averaged(_e , _l , Tref)
    end
end

function _pyrometer_band_averaged(λ1::NTuple{N}, λ2::NTuple{N}  , ϵ_func::AbstractSpectralQuantity , Tref) where N
    ntuple(N) do i
        l , r = λ1[i] , λ2[i]
        (emin , emax) = ϵ_func(l) , ϵ_func(r)
        ϵ_baseline = (emin + emax) / 2 # the idea is to exclude the average value 
        irem, _ = quadgk(λ -> (ϵ_func(λ , Tref) - ϵ_baseline) * Planck.ibb(λ, Tref), l, r)
        denom = Planck.band_power(Tref , λₗ = l, λᵣ = r )
        return (ϵ_baseline * denom + irem)/denom
    end
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
    set_emissivity!(p::AbstractPyrometer,  ϵ::AbstractContinuousOrDiscreteQuantity , Tref::Number)

Setting the emissivity from external data `ϵ`

# Arguments
`λ` - wavelength , 
`ϵ` - spectral emissivity 
"""
function set_emissivity!(p::AbstractPyrometer,  ϵ::AbstractContinuousOrDiscreteQuantity , Tref::Number)
    e_int = integral_emissivity(p ,  λ , ϵ , Tref)
    set_emissivity!(p , e_int)
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
    #Base.show(io::IO , p::MultiWavelengthPyrometer{}) 
    shorthand(p)=  "$(p.type) : $(p.λ), μm"
    include("custom_integration_and_interpolation_funcs.jl")
end