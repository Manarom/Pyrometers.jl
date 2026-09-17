
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
        convert_temperature,
        integral_emissivity,
        DefaultPyrometersTypes,
        fit_ϵ! , fit_ϵ , 
        Pyrometer , RatioPyrometer , 
        TabularQuantity , AnalyticalSpectralQuantity ,
        IsothermalSpectralQuantity , GenericDifferentiableSpectralQuantity , 
        stray_radiation_corrected_temperature , PlanckEmitter
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
    include("pyrometers_types.jl")
    abstract type AbstractDiscreteQuantity{LT , ET} end
    """
    subrange_view(λ1::Number , λ2::Number , e::AbstractDiscreteQuantity)

The view of subranges of discrete data within the `λ1...λ2` spectral range 

Must return the `Tuple{<:AbstractVector , <:AbstractVector} : (l  , e)`, where `l` and `e` are wavelength and quantity views 

"""
subrange_view(λ1::Number , λ2::Number , e::AbstractDiscreteQuantity) = error(" Undefined ")
    """
    subrange_view(λ1::NTuple{2} , λ2::NTuple{2} , e::AbstractDiscreteQuantity)

The view of subrange of discrete data within the `λ1[1]...λ1[2]` and `λ2[1]...λ2[2]` spectral range 

Must return the `NTuple{2 , D} where Tuple{<:AbstractVector , <:AbstractVector} : ((l1 , l2)  , (e1 , e2))` ,
 where `l1 - e1 `and `l2 - e2` are wavekengths and values views within the first and the second band 
"""
subrange_view(λ1::NTuple{2} , λ2::NTuple{2} , e::AbstractDiscreteQuantity) = error(" Undefined ")
    """
    get_single_wavelength_value(l::Number, _ , e::AbstractDiscreteQuantity ) -> Number
    get_single_wavelength_value(l::NTuple{N}, _ , e::AbstractDiscreteQuantity ) -> NTuple{N}

Functions to extract single interpolated values from disctrete data 
"""
function get_single_wavelength_value end

(tbq::AbstractDiscreteQuantity)(λ::Number) = get_single_wavelength_value(λ , nothing ,  tbq)
(tbq::AbstractDiscreteQuantity)(λ::Number , ::Number) = get_single_wavelength_value(λ , nothing ,  tbq)
integrate(l1::Number , l2::Number , i::AbstractDiscreteQuantity) = begin 
        (_l , _i) = subrange_view(l1 , l2 , i)
        return _simpson(_l , _i)
end
integrate(l1 , l2 , _::Number , i::AbstractDiscreteQuantity) = integrate(l1,l2 ,i)
"""
    TabularQuantity{LT <: AbstractVector , ET <: AbstractVector} <: AbstractDiscreteQuantity{LT , ET }

Type wrapper around discrete quantity with two columns `λ` and `i`

"""
    struct TabularQuantity{LT <: AbstractVector , ET <: AbstractVector} <: AbstractDiscreteQuantity{LT , ET }
        λ::LT
        i::ET
        TabularQuantity(l::LT , e::ET) where {LT <: AbstractVector , ET <: AbstractVector} = begin 
            @assert issorted(l) "First argument must be sorted in ascending order"
            @assert length(l) == length(e) "Two vectors must be of the same length"
            return new{LT , ET}(l , e)
        end
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
        return ( (l1 , l2) , (e1 , e2))
    end
    get_single_wavelength_value(l::Number , _ , e::TabularQuantity ) = _local_interpolate(l , e.λ , e.i)
    get_single_wavelength_value(l::NTuple{N} , _ , e::TabularQuantity ) where N = ntuple(N) do i 
        _local_interpolate(l[i] , e.λ , e.i)
    end 

    
    # continuous spectral quantities 
    const IsothermalSpectralQuantity = Planck.IsothermalSpectralQuantity 
    const AnalyticalSpectralQuantity = Planck.AnalyticalSpectralQuantity
    const AbstractSpectralQuantity =  Planck.AbstractSpectralQuantity
    get_single_wavelength_value(l::Number, t::Number , e::AbstractSpectralQuantity ) = e(l , t)
    get_single_wavelength_value(l::NTuple{N}, t::Number , e::AbstractSpectralQuantity) where N =ntuple(N) do i 
        e(l[i] , t)
    end



    Planck.eval_Dₜ(a::T , _ , _) where T <: Number = (a , zero(T) , zero(T))
    Planck.eval_Dₜ(a::T , _ ) where T <: Number = (a , zero(T) , zero(T))

    const ASQ = Planck.AbstractSpectralQuantity
    """
        Wrapper around two AbstractSpectralQuantity or Number product  
    """
    struct SpectralQuantitiesProduct{SQ1 , SQ2} <: ASQ
        e1::SQ1
        e2::SQ2
        function SpectralQuantitiesProduct(e1::SP1 , e2::SP2) where {SP1 <:Union{ASQ , Number} ,
                 SP2 <: ASQ}
                 return new{SP1 , SP2}(e1 , e2)
        end   
 
    end    
    SpectralQuantitiesProduct(e1::SP1 , e2::SP2) where {SP1 <: ASQ , SP2 <: Number} = SpectralQuantitiesProduct(e2 , e1)


    (qp::SpectralQuantitiesProduct{<:ASQ , <:ASQ})(λ , t) = qp.e1(λ , t) * qp.e2(λ , t) 
    (qp::SpectralQuantitiesProduct{<:Number , <:ASQ})(λ , t) = qp.e1 * qp.e2(λ , t) 
    (qp::SpectralQuantitiesProduct{<:Number , <:Number})(λ , t) = qp.e1 * qp.e2

    function Planck.eval_Dₜ(qp::SpectralQuantitiesProduct , l , t)
        (e1 , de1 , dde1) = Planck.eval_Dₜ(qp.e1 , l , t)
        (e2 , de2 , dde2) = Planck.eval_Dₜ(qp.e2 , l , t)
        return (e1 * e2 , 
                de1 * e2 + de2 * e1 , 
                dde1 * e2 + 2 * de1 * de2 + e1 * dde2)
    end 
  
    Base.:*(sq1::AbstractSpectralQuantity , sq2::AbstractSpectralQuantity) = SpectralQuantitiesProduct(sq1 , sq2)
    Base.:*(sq1::Number , sq2::AbstractSpectralQuantity) = SpectralQuantitiesProduct(sq1 , sq2)
    Base.:*(sq1::AbstractSpectralQuantity , sq2::Number) = SpectralQuantitiesProduct(sq1 , sq2)

    struct SpectralQuantitiesRatio{SQ1 , SQ2} <: AbstractSpectralQuantity
        e1::SQ1
        e2::SQ2
        function SpectralQuantitiesRatio(e1::SP1 , e2::SP2) where {SP1<: Union{AbstractSpectralQuantity , Number} ,
                 SP2 <: Union{AbstractSpectralQuantity , Number} }
                 return new{SP1 , SP2}(e1 , e2)
        end   
    end    

    
    (qp::SpectralQuantitiesRatio)(λ , t) = qp.e1(λ , t) / qp.e2(λ , t) 
    (qp::SpectralQuantitiesRatio{<:Number})(λ , t) = qp.e1 / qp.e2(λ , t) 
    (qp::SpectralQuantitiesRatio{<:ASQ , <:Number})(λ , t) = qp.e1(λ , t)  / qp.e2
    (qp::SpectralQuantitiesRatio{<:Number , <:Number})(λ , t) = qp.e1 / qp.e2

    function Planck.eval_Dₜ(qp::SpectralQuantitiesRatio , l , t)
        (e1 , de1 , dde1) = Planck.eval_Dₜ(qp.e1 , l , t)
        (e2 , de2 , dde2) = Planck.eval_Dₜ(qp.e2 , l , t)
        return (
                e1 / e2 , 
                Planck._spectral_ratio_first_derivative(e1 , de1 , e2 , de2), 
                Planck._spectral_ratio_second_derivative(e1 , de1 , dde1 , e2 , de2 , dde2)
        )
    end     
    Base.:/(sq1::Union{AbstractSpectralQuantity , Number} , sq2::AbstractSpectralQuantity ) = SpectralQuantitiesRatio(sq1 , sq2)
    Base.:/(sq1::AbstractSpectralQuantity , sq2::Union{AbstractSpectralQuantity , Number}  ) = SpectralQuantitiesRatio(sq1 , sq2)

    struct SpectralQuantitiesSum{S1 , S2} <: AbstractSpectralQuantity
        e1::S1
        e2::S2
        function SpectralQuantitiesSum(e1::SP1 , e2::SP2) where {SP1<: Union{AbstractSpectralQuantity , Number} ,
                 SP2 <: Union{AbstractSpectralQuantity , Number}}
                 return new{SP1 , SP2}(e1 , e2)
        end   
    end    
    (sqs::SpectralQuantitiesSum)(l , t) = sqs.e1(l , t) + sqs.e2(l , t)
    (sqs::SpectralQuantitiesSum{S1})(l , t) where {S1<: Number} = sqs.e1 + sqs.e2(l , t)
    Planck.eval_Dₜ(sqs::SpectralQuantitiesSum , l , t) = begin 
        (e1 , de1 , dde1) = Planck.eval_Dₜ(sqs.e1 , l , t)
        (e2 , de2 , dde2) = Planck.eval_Dₜ(sqs.e2 , l , t)
        return (e1 + e2 , de1 + de2 , dde1 + dde2)
    end
    Planck.eval_Dₜ(sqs::SpectralQuantitiesSum{S1} , l , t) where {S1<:Number} = begin 
        (e1 , de1 , dde1) = Planck.eval_Dₜ(sqs.e1 , l , t)
        (e2 , de2 , dde2) = Planck.eval_Dₜ(sqs.e2 , l , t)
        return (e1 + e2 , de1 + de2 , dde1 + dde2)
    end
    Base.:+(asq1::ASQ , asq2::ASQ) = SpectralQuantitiesSum(asq1 , asq2)
    Base.:+(asq1::Number , asq2::ASQ) = SpectralQuantitiesSum(asq1 , asq2)
    Base.:+(asq1::ASQ , asq2::Number) = SpectralQuantitiesSum(asq2 , asq1)
    Base.:-(asq1::ASQ , asq2::Number) = SpectralQuantitiesSum(-asq2 , asq1)
    Base.:-(asq1::T , asq2::ASQ) where T <: Number = SpectralQuantitiesSum(asq1 , -one(T) * asq2)
    Base.:-(asq1::ASQ, asq2::ASQ)  = SpectralQuantitiesSum(asq1 , (-1.0) * asq2)


    struct PlanckEmitter <: AbstractSpectralQuantity   end
    (::PlanckEmitter)(l , t) = Planck.ibb(l , t)
    Planck.eval_Dₜ(::PlanckEmitter , l , t) = Planck.Dₜibb(l , t)
    """
    fix_temperature(q::AbstractSpectralQuantity , fixed_temperature::Number) -> ::IsothermalSpectralQuantity

Fixes `AbstractSpectralQuantity` temperature converting it to `IsothermalSpectralQuantity`
"""
fix_temperature(q::AbstractSpectralQuantity , fixed_temperature::Number) = IsothermalSpectralQuantity(Base.Fix2(q , fixed_temperature))


    (p::Planck.IsothermalSpectralQuantity)(λ) = p(λ , nothing)
   
    integrate(l1::Number , l2::Number , t::Number ,  e::E; segbuf=nothing, kwargs...) where E <: AbstractSpectralQuantity =quadgk(Base.Fix2(e , t) , l1 , l2 ; segbuf=segbuf, kwargs...)[1]
    integrate(l1::Number , l2::Number  ,  e::E; segbuf=nothing, kwargs...) where E <: IsothermalSpectralQuantity = quadgk(e , l1 , l2 ; segbuf=segbuf , kwargs...)[1] 


    struct SpectralQuantityIntegrator{SQ , L}
        sq::SQ
        l1::L 
        l2::L
        SpectralQuantityIntegrator(sq::SQ , l1::L , l2::L) where {L <: Number , SQ <: AbstractSpectralQuantity} = new{SQ , L}(sq , l1 , l2)
    
    end
    #function SpectralQuantityIntegrator(p::SpectralBandPyrometer , sq::AbstractSpectralQuantity) 
    #    SpectralQuantityIntegrator(sq , p.λ[1] , p.λ[2])
    #end 
    (sqi::SpectralQuantityIntegrator)(t::Number; kwargs...) = integrate(sqi.l1 , sqi.l2 , t , sqi.sq; kwargs...)
    function Planck.eval_Dₜ(sqi::SpectralQuantityIntegrator , t::Number) 
        f(l) = SVector(Planck.eval_Dₜ(sqi.sq , l , t))
        return Tuple(first(quadgk(f , sqi.l1 , sqi.l2)))
    end
    struct SpectralQuantityIntegratorContext{SQI , T}
        sqi::SQI
        i::T 
    end
    """
    fit_integral(a::AbstractSpectralQuantity , l1, l2 , measured)

Fits the temperature of any spectral quantity to a specified value  , 
"""
function fit_integral(a::AbstractSpectralQuantity , l1, l2 , measured) 
    sqi = SpectralQuantityIntegrator(a  ,l1 , l2)
    return fit_spectral_quantity_integrator(sqi , measured)
end
function fit_spectral_quantity_integrator(sqi::SpectralQuantityIntegrator , imeasured::Number; T_starting=600.0)
        ctx = SpectralQuantityIntegratorContext(sqi, imeasured)
        return Roots.find_zero(ctx , T_starting ,  Roots.Halley())
    end
    function (sqic::SpectralQuantityIntegratorContext)(t) 
        _to_halley(Planck.eval_Dₜ(sqic.sqi , t) , sqic.i)
    end    
    Planck.∫ₗ(a::AbstractSpectralQuantity , l1 , l2) = SpectralQuantityIntegrator(a , l1 ,l2)
    const AbstractContinuousOrDiscreteQuantity = Union{AbstractSpectralQuantity , AbstractDiscreteQuantity}
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


include("measure_func.jl")

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
    """
    convert_temperature(p::AbstractPyrometer , Tmeasured  , ϵ_new)

Converts temperature `Tmeasured` measured using pyrometer `p` with it specified emissivity 
to a new temperature measured with `ϵ_new` , the type of `ϵ_new` depends on the type of pyrometer 
if `ϵ_new` is a `Number` than if p is `RatioPyrometer` it assumes `e_new` is `e_slope`, if 
`e_new` is `NTuple{2 , Number}` it modifies both emissivities at two wavelength
"""
convert_temperature(p::AbstractPyrometer , Tmeasured  ,  ϵ_new::Union{Number , AbstractContinuousOrDiscreteQuantity}) = measure(p , signal(p , Tmeasured) , ϵ_new)

convert_temperature(p::AbstractPyrometer , Tmeasured  ,  ϵ_previous::Union{Number , AbstractContinuousOrDiscreteQuantity}, 
                    ϵ_new::Union{Number , AbstractContinuousOrDiscreteQuantity}) = measure(p , signal(p , Tmeasured , ϵ_previous) , ϵ_new)


"""
    integrate(p::AbstractPyrometer ,  temperature:: Number , intensity::AbstractContinuousOrDiscreteQuantity) 
    integrate(p::SpectralBandPyrometer , 
                intensity::Union{AbstractDiscreteQuantity , IsothermalSpectralQuantity})
Converts the input spectral `intensity` given as a discrete set of values for wavelengths `λ`
to the pyrometer signal.

Returns the quantity, which is equal to the type of pyrometer signal , e.g. if pyrometer is 

- `SingleWavelengthPyrometer` returns the spectral intensity at the wavelength of pyrometer

- `SpectralBandPyrometer`  - total intensity within the spectral range the pyrometer

- `TwoWavelengthRatioPyrometer`  - tuple of intensities at two wavelengths 

- `TwoBandsRatioPyrometer`  - tuple of total intensities within two bands 


`λ` , `intensity`  (both at the same time) if the intensity is provided as a discrete set of points 

"""
@inline integrate(p::SpectralBandPyrometer , 
                intensity::Union{AbstractDiscreteQuantity , IsothermalSpectralQuantity}; segbuf=nothing, kwargs...) = integrate(p.λ[1] , p.λ[2] , intensity; segbuf=segbuf, kwargs...)

@inline integrate(p::SpectralBandPyrometer , t::Number , 
    intensity::AbstractSpectralQuantity; segbuf=nothing, kwargs...) = integrate(p.λ[1] , p.λ[2] , t , intensity; segbuf=segbuf, kwargs...)


@inline function integrate(p::TwoBandsRatioPyrometer , t,
     intensity_function::Union{AbstractDiscreteQuantity , IsothermalSpectralQuantity}; segbuf=nothing, kwargs...) 
    band1, band2 = p.λ[1], p.λ[2]
    i1 = integrate(band1[1] , band1[2] , t ,  intensity_function; segbuf=segbuf, kwargs...)
    i2 = integrate(band2[1] , band2[2] , t ,  intensity_function; segbuf=segbuf, kwargs...)
    return (i1 , i2)
end
@inline  function integrate(p::TwoBandsRatioPyrometer ,
     intensity_function::Union{AbstractDiscreteQuantity , IsothermalSpectralQuantity}; segbuf=nothing,  kwargs...) 
    band1, band2 = p.λ[1], p.λ[2]
    i1 = integrate(band1[1] , band1[2] ,  intensity_function; segbuf=segbuf, rtol=rtol , kwargs...)
    i2 = integrate(band2[1] , band2[2] ,  intensity_function; segbuf=segbuf, rtol=rtol , kwargs...)
    return (i1 , i2)
end
# versions for single wavelength pyrometers 
integrate(p::SingleWavelengthPyrometer  , intensity::Union{AbstractDiscreteQuantity , IsothermalSpectralQuantity}; segbuf=nothing , kwargs...) = intensity(p.λ[1])
integrate(p::SingleWavelengthPyrometer, t , intensity_function::AbstractContinuousOrDiscreteQuantity; segbuf=nothing , kwargs...) = intensity_function(p.λ[1] , t )
integrate(p::TwoWavelengthRatioPyrometer, intensity_function; segbuf=nothing,kwargs...) = (intensity_function(p.λ[1]) , intensity_function(p.λ[2]))
integrate(p::TwoWavelengthRatioPyrometer , t , intensity::AbstractContinuousOrDiscreteQuantity; segbuf=nothing, kwargs...) = begin 
    return (
            intensity(p.λ[1] , t) , 
            intensity(p.λ[2] , t)
            )
end
                    
abstract type AbstractRadiationGeometry end
""" 
    Type stores the geometry parameter `ξ= F * A₁/A₂` , where `F` is view factor  , 
    `A₁` and `A₂` are surface under measurement area and external heater area 
"""
struct ViewFactorGeometry{F <: Number} <: AbstractRadiationGeometry
         ξ::F # F_12 *A1/A2
         ViewFactorGeometry(ξ::F) where F <: Number = new{F}(ξ)
end
struct EnclosureGeometry <: AbstractRadiationGeometry end
ViewFactorGeometry(A1 , A2 , F12) = ViewFactorGeometry(F12 * A1/A2)
#EnclosureGeometry() = ViewFactorGeometry(0.0)
#ParallelGeometry() = ViewFactorGeometry(1.0)  

 """
    effective_emissivity(geometry::AbstractRadiationGeometry, eo::Number, es::Number)

returns effective emissivity 
"""
@inline function effective_emissivity(geometry::AbstractRadiationGeometry, eo::Number, es::Number)
    ξ = geometry.ξ
    denom = one(ξ) + eo * (one(es)/es - one(es)) * ξ
    ϵ_priv = eo / denom
    return ϵ_priv 
end

struct SpectralReflectivity{ S} <: AbstractSpectralQuantity
    ϵ::S
end

@inline function (r::SpectralReflectivity{S})(λ, t) where S <: AbstractContinuousOrDiscreteQuantity
    e_val = r.ϵ(λ , t)
    return one(e_val) - e_val
end

(r::SpectralReflectivity{S})(λ) where S <: Union{IsothermalSpectralQuantity , TabularQuantity , Number} = r(λ , nothing)

@inline (r::SpectralReflectivity{S})(λ, t::T) where {S <: Number , T}  = one(T) - r.ϵ

function Planck.eval_Dₜ(r::SpectralReflectivity{S} , λ , T) where S <: Union{AbstractContinuousOrDiscreteQuantity}
    e , de , dde  = Planck.eval_Dₜ(r.ϵ , λ , T)
    return (one(e) - e , -de ,-dde )
end
Planck.eval_Dₜ(r::SpectralReflectivity{<:Number} , λ , t::DT) where DT = (r(λ , t) , zero(DT) , zero(DT))

struct EffectiveEmissivityQuantity{ O, S, Tsrc <: Number , G <: AbstractRadiationGeometry} <: AbstractSpectralQuantity
    ϵ_surf::O
    ϵ_src::S
    Tsource::Tsrc 
    geom::G
end

@inline function (q::EffectiveEmissivityQuantity{O , S})(λ, t) where {O <: AbstractContinuousOrDiscreteQuantity , S <: AbstractContinuousOrDiscreteQuantity}
    eo = q.ϵ_surf(λ, t)          # measurement surface spectral emissivity  
    es = q.ϵ_src(λ, q.Tsource)   # source spectral emissivity 
    return effective_emissivity(q.geom, eo, es)
end
@inline function (q::EffectiveEmissivityQuantity{O , S})(λ, t) where {O <: AbstractContinuousOrDiscreteQuantity , S <: Number}
    eo = q.ϵ_surf(λ, t)          # measurement surface spectral emissivity  
    es = q.ϵ_src                 # source spectral emissivity 
    return effective_emissivity(q.geom, eo, es)
end
@inline function (q::EffectiveEmissivityQuantity{O , S})(_, _) where {O <: Number , S <: Number}
    eo = q.ϵ_surf        # measurement surface spectral emissivity  
    es = q.ϵ_src                 # source spectral emissivity 
    return effective_emissivity(q.geom, eo, es)
end
@inline function (q::EffectiveEmissivityQuantity{O , S})(λ, t) where {O <: Number , S <: AbstractContinuousOrDiscreteQuantity}
    eo = q.ϵ_surf        # measurement surface spectral emissivity  
    es = q.ϵ_src(λ, q.Tsource)   # source spectral emissivity 
    return effective_emissivity(q.geom, eo, es)
end
"""
    combine_effective_epsilon_derivatives(geom::ViewFactorGeometry, eo_tpl, es)

Converts the derivatives of effective emissivity according to the chain rule.
`eo_tpl` is a tuple of (value, first derivative, second derivative) of the target surface emissivity.
`es` is fixed and does not depend on temperature (but can depend on wavelength).
"""
@inline function combine_effective_epsilon_derivatives(geom::ViewFactorGeometry, eo_d3 ,  es)
    eo, d_eo, dd_eo = eo_d3
    ξ = geom.ξ
    k  = ξ * (one(es)/es  - 1) # this coefficient does not depend on temperature 
    k_k = (1 + eo * k )
    denom = one(eo)/k_k
    e_eff = eo * denom 
    denom_sqr = denom ^2 
    # e_eff =  eo / (1 + eo * k )
    # d_e_eff = (d_eo * ( 1 + eo * k) - eo * d_eo * k)/(1 + eo * k)² = d_eo /(1 + eo * k)²  
    # u = eo * es , effective emissivity first derivative 
    d_e_eff = d_eo * denom_sqr
    dd_e_eff = (dd_eo  - 2 * d_eo * d_eo * k * denom) * denom_sqr    
    return (e_eff , d_e_eff , dd_e_eff)
end
function Planck.eval_Dₜ(e_eff::EffectiveEmissivityQuantity{<:AbstractContinuousOrDiscreteQuantity , <:AbstractContinuousOrDiscreteQuantity} , λ , T)
    eo_d3  = Planck.eval_Dₜ(e_eff.ϵ_surf , λ , T)
    return combine_effective_epsilon_derivatives(e_eff.geom , eo_d3 , e_eff.ϵ_src(λ , e_eff.Tsource))
end
function Planck.eval_Dₜ(e_eff::EffectiveEmissivityQuantity{<:AbstractContinuousOrDiscreteQuantity , <:Number} , λ , T)
    eo_d3  = Planck.eval_Dₜ(e_eff.ϵ_surf , λ , T)
    return combine_effective_epsilon_derivatives(e_eff.geom , eo_d3 , e_eff.ϵ_src)
end
function Planck.eval_Dₜ(e_eff::EffectiveEmissivityQuantity{ST, <:AbstractContinuousOrDiscreteQuantity} , λ , T) where ST <: Number
    eo_d3  = (e_eff.ϵ_surf , zero(ST) , zero(ST))
    return combine_effective_epsilon_derivatives(e_eff.geom , eo_d3 , e_eff.ϵ_src(λ , e_eff.Tsource))
end
function Planck.eval_Dₜ(e_eff::EffectiveEmissivityQuantity{<: Number, <:Number} , λ , T::D) where D <: Number  
    return (e_eff(λ , T) , zero(D) , zero(D))
end

#basic types functions 
"""
    stray_radiation_corrected_temperature(p::Pyrometer, Tmeasured::Number, source_intensity::Number)

Correct the measured temperature for a single-channel brightness pyrometer when the external 
background radiation (`source_intensity`) is provided directly as an integrated signal.

The baseline surface emissivity is automatically retrieved from the pyrometer instance. 
The system geometry is assumed to be an enclosure (`EnclosureGeometry`), meaning multiple mutual 
reflections are neglected, yielding a simple effective reflectivity of `r_eff = 1 - ϵ_surf`.
"""
function stray_radiation_corrected_temperature(p::Pyrometer, Tmeasured::Number, 
                                                    ϵ_surf::Number,
                                                    source_intensity::Number) 

        measured_signal = signal(p, Tmeasured , ϵ_surf)         
        r_eff   = one(ϵ_surf) - ϵ_surf
        reflected_signal = r_eff * source_intensity
        pure_signal = _extract_signals(p , measured_signal , reflected_signal)
        return p(pure_signal , ϵ_surf)
    end

"""
    stray_radiation_corrected_temperature(p::Pyrometer, Tmeasured::Number, 
                                                    ϵ_surf::AbstractSpectralQuantity,
                                                    source_intensity::Union{Number, IsothermalSpectralQuantity})

Surface emissivity can be temperature-dependent , external radiation must be provided as an isothermal quantity 
"""
function stray_radiation_corrected_temperature(p::Pyrometer, Tmeasured::Number, 
                                                    ϵ_surf::AbstractSpectralQuantity,
                                                    source_intensity::Union{Number, IsothermalSpectralQuantity}) 

        measured_signal = signal(p, Tmeasured , ϵ_surf)         
        r_eff   = SpectralReflectivity(ϵ_surf)
        reflected_signal = r_eff * source_intensity
        pure_signal = fix_temperature(_extract_signals(p , measured_signal , reflected_signal) , Tmeasured)
        return p(pure_signal , ϵ_surf) #
    end


"""
    stray_radiation_corrected_temperature(p::RatioPyrometer, Tmeasured::Number, 
                                                    source_intensity::NTuple{2 , T}) where T <: Number

Correct the measured temperature for a dual-channel ratio pyrometer when the external 
background radiation signals for both channels are directly provided as a tuple (`source_intensity`).

The baseline surface emissivities for both channels are automatically retrieved from the pyrometer instance. 
The system geometry is assumed to be an enclosure (`EnclosureGeometry`), meaning multiple mutual 
reflections are neglected, yielding independent channel reflectivities `r_eff = 1 - e_surf`. 
Returns the true temperature resolved from the corrected signal ratio.
"""
function stray_radiation_corrected_temperature(p::RatioPyrometer, Tmeasured::Number, e_surf::NTuple{2},
                                                    source_intensity::NTuple{2 , T}) where T <: Number

        e_surf1 , e_surf2 = e_surf
        si1 , si2 = source_intensity[1] , source_intensity[2]
        m_s1 , m_s2 = signal(p, Tmeasured , e_surf)         
        r_eff1 , r_eff2   = ( one(e_surf1) - e_surf1 , one(e_surf2) - e_surf2)
        r_s1 , r_s2 = r_eff1 * si1 , r_eff2 * si2 
        p_s1 , p_s2 = m_s1 - r_s1 , m_s2 - r_s2
        return p(p_s1/p_s2)
    end
"""
    stray_radiation_corrected_temperature(p::AbstractPyrometer, 
                                        Tmeasured::Number, 
                                        source_intensity
                                        ) where T<:Number

Takes emissivity from the pyrometer 
"""
stray_radiation_corrected_temperature(p::AbstractPyrometer, 
                                        Tmeasured::Number, 
                                        source_intensity
                                        )  = stray_radiation_corrected_temperature(p, Tmeasured,  
                                                get_emissivity(p) ,
                                                source_intensity
                                            )

"""
    stray_radiation_corrected_temperature(p::AbstractPyrometer, 
                                        Tmeasured::Number, 
                                        ϵ_surf::Union{T, NTuple{2, T}} ,
                                        source_intensity::AbstractContinuousOrDiscreteQuantity 
                                        ) where T<:Number

when the surface emissivity is a `Number` external radiation can be integrated 
"""
function stray_radiation_corrected_temperature(p::AbstractPyrometer, 
                                        Tmeasured::Number, 
                                        ϵ_surf::Union{T, NTuple{2, T}} ,
                                        source_intensity::Union{AbstractDiscreteQuantity , IsothermalSpectralQuantity};
                                        segbuf=nothing, rtol=sqrt(eps(Float64))
                                        ) where T<:Number 
                                        
    return stray_radiation_corrected_temperature(p, Tmeasured,  
                                                ϵ_surf ,
                                                integrate(p , Tmeasured , source_intensity , segbuf=segbuf, rtol=rtol)
                                        )
end
"""
    stray_radiation_corrected_temperature(p::AbstractPyrometer, 
        Tmeasured::Number, 
        ϵ_surf::Union{T , NTuple{2,T}},
        incident_radiation_function::AbstractSpectralQuantity ,  
        radiation_temperature::Number ; 
        segbuf=nothing, rtol=sqrt(eps(Float64)) ) where {T <: Number}

Correct the measured temperature when the external background source emits a known temperature-dependent 
intensity profile (which is not required to be an ideal blackbody spectrum).

"""
function stray_radiation_corrected_temperature(p::AbstractPyrometer, 
                                                Tmeasured::Number, 
                                                ϵ_surf::Union{T , NTuple{2,T}},
                                                incident_radiation_function::AbstractSpectralQuantity ,  
                                                radiation_temperature::Number ; 
                                                segbuf=nothing, rtol=sqrt(eps(Float64)) ) where {T <: Number}
        
        return stray_radiation_corrected_temperature(p , Tmeasured ,
                    ϵ_surf, 
                    integrate(p  , radiation_temperature ,  
                                incident_radiation_function , 
                                segbuf=segbuf, rtol=rtol)
        )
    end  




@inline function _extract_signals(::Pyrometer, measured, reflected)
    return measured - reflected
end

@inline function _extract_signals(::RatioPyrometer, measured::NTuple{2}, reflected::NTuple{2})
    m1, m2 = measured
    r1, r2 = reflected
    return (m1 - r1) / (m2 - r2)
end
    

  


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
    
    shorthand(p)=  "$(p.type) : $(p.λ), μm"
    include("custom_integration_and_interpolation_funcs.jl")
end