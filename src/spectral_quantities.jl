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
integrate(l1 , l2 , _::Number , i::AbstractDiscreteQuantity) = integrate(l1 , l2 , i)
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
    Base.:/(sq1::Number , sq2::ASQ ) = SpectralQuantitiesRatio(sq1 , sq2)
    Base.:/(sq1::ASQ , sq2:: Number  ) = SpectralQuantitiesRatio(sq1 , sq2)
    Base.:/(sq1::ASQ , sq2::ASQ ) = SpectralQuantitiesRatio(sq1 , sq2)

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
    """
        Common type for all continuous or discrete quantities 
    """
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