
module RadiationPyrometers
       #Optimization,
            #OptimizationOptimJL,
    using   LinearAlgebra,
            StaticArrays,
            OrderedCollections, 
            Roots
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
                           Ref(1.0)) 
            else
                 error("Unknown pyrometer type")
            end
        end
        Pyrometer(λ::NTuple{N , T}; type::Symbol=:def ,  ϵ::Number=1.0 ) where {N,T} = begin 
            return new{N , T}(type,SVector{N}(λ),Ref(ϵ))
        end
        Pyrometer(λ::Union{AbstractVector{T} , T} ; type::Symbol=:def ,  ϵ::Number=1.0) where T <: Number = begin
            #@assert 0  < ϵ <= 1.0 "Emissivity should be within the (0..1] interval"
            if λ isa AbstractVector
                N = length(λ)
                @assert N == 1 || N == 2 "λ should be a vector of two  Floats or a single Float number"
            else
                N = 1
            end
            new{N , T}(type,SVector{N}(λ),Ref(ϵ))
        end
    end
    """
    wlength(::Pyrometer{N}) where N

Returns the number of wavelengths
"""
wlength(::Pyrometer{N}) where N = N
is_narrow_band(::AbstractPyrometer)  = false
    """
    is_narrow_band(p::Pyrometer)

True if pyrometer `p` is a narrow-band pyrometer (works on a fixed wavelengh region)
"""
is_narrow_band(::Pyrometer{2})  = true
is_narrow_band(::RatioPyrometer{2 , T , DT}) where {T , DT <: Tuple}  = true
"""
    is_fixed_wavelength(::Pyrometer{N}) where N

True if Pyrometer is single wavelength
"""
is_fixed_wavelength(::Pyrometer{1}) = true
is_fixed_wavelength(::RatioPyrometer{2 , T , DT}) where {T , DT <: Number} = true
is_fixed_wavelength(::AbstractPyrometer) = false
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
        measure_with_env(p::Pyrometer{2}, i_total::Number, Tenv::Number; T_starting=600.0)

    
    """
    function measure_with_env(p::Pyrometer, Tmeasured::T, Tenv::Number , ϵ_env::Number = 1.0) where {T}
        ϵ = _get_epsilon_equivalent(p)
        measured_signal = signal(p , Tmeasured)
        reflected_signal = ϵ_env * (one(T) - ϵ) * signal(p , Tenv) / ϵ
        return p(measured_signal - reflected_signal)
    end    
    """
    signal(p::Pyrometer , Tmeasured)

Returns the signal value which will give the temperature `Tmeasured`
"""
signal(p::Pyrometer{1} , Tmeasured::Number)  = p.ϵ[] * Planck.ibb(p.λ[] , Tmeasured)
signal(p::Pyrometer{2} , Tmeasured::Number) = p.ϵ[] * Planck.band_power(Tmeasured , λₗ=p.λ[1] , λᵣ=p.λ[2])
signal(p::RatioPyrometer{2 , T , DT} , Tmeasured::Number)  where { T , DT <: Number} = Planck.spectral_ratio(p.λ[1] , p.λ[2] , Tmeasured , e_slope = e_slope(p))
signal(p::RatioPyrometer{2 , T , DT} , Tmeasured::Number)  where { T , DT <: Tuple} = Planck.spectral_band_ratio(p.λ[1] , p.λ[2] , Tmeasured , e_slope = e_slope(p))

    """
    fit_ϵ(p::AbstractPyrometer , Tmeasured::Number , Treal::Number)

Finds the emissivity or e_slope 
Input:
p - pyrometer object
Treal - real temperature of the surface, Kelvins
Tmeasured - temperature measured by the pyrometer, Kelvins
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
            if is_narrow_band(pj)
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
    set_emissivity!(p::Pyrometer,em_value::Float64)

Setter for spectral emissivity
"""
set_emissivity!(p::Pyrometer , em_value::Number) = (p.ϵ[] = em_value)
set_emissivity!(p::RatioPyrometer , em_value::Number) = begin 
    p.ϵ1[]  = em_value * p.ϵ2[] 
end
set_emissivity!(p::RatioPyrometer{N , T} , em_value::NTuple{2 , D}) where {D <: Number , N , T} = begin 
    p.ϵ1[]  = T(em_value[1])
    p.ϵ2[]  = T(em_value[2]) 
end

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
            if is_narrow_band(p[i]) 
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
function fit_ϵ_wavelength!(p::Pyrometer,Tmeasured::Float64,Treal::Float64) # this is the same as fit_ϵ! with the exception that 
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
    Base.show(io::IO, p::Pyrometer{1}) = print(io, "$(p.type) -type: Fixed-wavelength pyrometer:λ = $(p.λ[1]) μm,ϵ = $(p.ϵ[])")
    Base.show(io::IO, p::Pyrometer{2}) = print(io, "$(p.type) -type: Narrow-band pyrometer:λ ∈ $(p.λ[1]) ... $(p.λ[2]) μm,ϵ = $(p.ϵ[])")
    Base.show(io::IO, p::RatioPyrometer{2 , T , DT}) where {T , DT <: Number} = print(io, "$(p.type) -type: Fixed-wavelength's spectral ratio pyrometer:λ₁= $(p.λ[1]) , λ₂ = $(p.λ[2]) μm, ϵ₁ = $(p.ϵ1[]) , ϵ₂ = $(p.ϵ2[]) , e_slope = $(e_slope(p))")
    Base.show(io::IO, p::RatioPyrometer{2 , T , DT}) where {T , DT <: Tuple} = print(io, "$(p.type) -type: Narrow-band spectral ratio pyrometer:λ₁= $(p.λ[1]) , λ₂ = $(p.λ[2]) μm, ϵ₁ = $(p.ϵ1[]) , ϵ₂ = $(p.ϵ2[]) , e_slope = $(e_slope(p))")
end