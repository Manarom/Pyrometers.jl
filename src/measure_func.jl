
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
`T_starting`  - temperature hint

"""
measure(p::AbstractPyrometer , i::Number  , ϵ::Union{Number , NTuple{2,<:Number}}; 
                T_starting::Number=1000.0  , 
                segbuf=nothing, kwargs...) = Roots.find_zero(t -> _Dₜpyro(p.λ , i , t , ϵ) , T_starting ,  Roots.Halley() ; kwargs...) 

measure(p::RatioPyrometer , i::NTuple{2,D} ,  ϵ::NTuple{2 , <:Number} ; kwargs...) where {D <: Number} = measure(p , i[1]/i[2] , ϵ ; kwargs...)


"""
    measure(p::AbstractPyrometer , i::Union{IsothermalSpectralQuantity , AbstractDiscreteQuantity}; 
                    T_starting::Number = 1000.0  , 
                    segbuf=nothing, kwargs...)

Input radiation is provided as temperature independent spectral quantity  , the emissivity is taken from the 
pyrometer's default (thus we know it is both wavelength and temperature independent,probably integrated)
"""
measure(p::AbstractPyrometer , i::Union{IsothermalSpectralQuantity , AbstractDiscreteQuantity}; 
                    T_starting::Number = 1000.0  , 
                    segbuf=nothing, kwargs...) = measure(p , integrate(p , i; segbuf=segbuf, kwargs...) , T_starting = T_starting)

"""
    measure(p::AbstractPyrometer , i::AbstractContinuousOrDiscreteQuantity,
                        radiation_temperature::Number; segbuf=nothing, kwargs...)

Input radiation intensity is provided as a temperature-dependent spectral quantity with known temperature 
(it can be some sort of parameter, not temperature) 
"""
measure(p::AbstractPyrometer , i::AbstractContinuousOrDiscreteQuantity,
                        radiation_temperature::Number; segbuf=nothing, kwargs...) = measure(p , integrate(p , radiation_temperature ,  
                                                                                    i ; segbuf=segbuf , kwargs...) ; kwargs...)

function measure(p::AbstractPyrometer , i::Union{Number , NTuple{2,<:Number}} ; kwargs...)
        _ϵ = _get_epsilon_equivalent(p)
        return measure(p , i , _ϵ; kwargs...)
    end

Dₜpyro(p::AbstractPyrometer , i , t) = _Dₜpyro(p.λ , i , t , _get_epsilon_equivalent(p))
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
    measure(p::AbstractPyrometer  , 
                    imeasured::Number , 
                    ϵ::AbstractContinuousOrDiscreteQuantity; 
                    T_starting::Number = 600.0)

    Emissivity of the surface is provided as spectral quantity but measured intencity is a number
    (somebody else integrated it)
    
# Arguments
- `p`: AbstractPyrometer object 
- `imeasured` : measured signal 
- `ϵ` : surface emissivity, can be provided as 
(optional)
- `T_starting` : starting temperature for nonlinear eqaution solver

"""
function measure(p::AbstractPyrometer  , 
                    imeasured::Number , 
                    ϵ::AbstractContinuousOrDiscreteQuantity; 
                    T_starting::Number = 600.0 ,  segbuf=nothing, kwargs...)

    ctx = SpectralQuantityPyrometricContext(p , ϵ , imeasured)
    return Roots.find_zero(ctx , T_starting ,  Roots.Halley() , kwargs...)
end
"""
    measure(p::RatioPyrometer  , 
                    imeasured::NTuple{2} , 
                    ϵ::AbstractContinuousOrDiscreteQuantity; 
                    T_starting::Number = 600.0 ,  segbuf=nothing, rtol=sqrt(eps(Float64)))

Evaluates the temperature from signal measured in two channels of two-color pyrometer 
"""
measure(p::RatioPyrometer  , 
                    imeasured::NTuple{2} , 
                    ϵ::AbstractContinuousOrDiscreteQuantity; 
                    T_starting::Number = 600.0 ,  segbuf=nothing, kwargs...) = measure(p , imeasured[1]/imeasured[2] , ϵ ; T_starting = T_starting ,   segbuf=segbuf, kwargs...)

"""
    measure(p::AbstractPyrometer , i::Union{IsothermalSpectralQuantity , AbstractDiscreteQuantity} , 
                    ϵ::AbstractContinuousOrDiscreteQuantity; 
                    T_starting::Number = 600.0 , segbuf=nothing, rtol=sqrt(eps(Float64)))


Measures temperature from external radiation `i` provided as a spectral quantity which does'n depend
on temperature, taking into account the surface emissivity `ϵ`
"""
measure(p::AbstractPyrometer , i::Union{IsothermalSpectralQuantity , AbstractDiscreteQuantity} , 
                    ϵ::AbstractContinuousOrDiscreteQuantity; 
                    T_starting::Number = 600.0 , segbuf=nothing , kwargs...) = measure(p , integrate(p , i ; segbuf=segbuf , kwargs...) , ϵ , T_starting = T_starting , kwargs...)
"""
    measure(p::AbstractPyrometer , i::Union{IsothermalSpectralQuantity , AbstractDiscreteQuantity} , 
                    ϵ::AbstractContinuousOrDiscreteQuantity; 
                    T_starting::Number = 600.0)

Measures temperature from external radiation `i` provided as a spectral quantity which depends on temperature 
`radiation_temperature`, taking into acoount the surface emissivity `ϵ` , which can depend on surface temperature
"""
measure(p::AbstractPyrometer , i::AbstractContinuousOrDiscreteQuantity ,
                    radiation_temperature::Number , 
                    ϵ::AbstractContinuousOrDiscreteQuantity; 
                    T_starting::Number = 600.0 , segbuf=nothing, rtol=sqrt(eps(Float64))) = measure(p , integrate(p , radiation_temperature , i; segbuf=segbuf, rtol=rtol) , ϵ ; T_starting = T_starting)

## functors 

struct SpectralQuantityPyrometricContext{L , E , F , P} 
    λ::L
    e::E
    i_measured::F
    p::P
    function SpectralQuantityPyrometricContext(p::P , quantity::E , i_measured::F) where { E , F , P<:AbstractPyrometer} 
        _λ, _e = _prepare_context_data(p, quantity)
        return new{typeof(_λ), typeof(_e), F, P}(_λ, _e, i_measured, p)
    end
end

struct GenericSpectralQuantityIntegrationContext{L , E , F , P} 
    λ::L
    e::E
    i_measured::F
    p::P
    function GenericSpectralQuantityIntegrationContext(p::P , quantity::E , i_measured::F) where { E , F , P<:AbstractPyrometer} 
        _λ, _e = _prepare_context_data(p, quantity)
        return new{typeof(_λ), typeof(_e), F, P}(_λ, _e, i_measured, p)
    end
end


@inline function _prepare_context_data(p::AbstractPyrometer, quantity::AbstractSpectralQuantity)
    return Tuple(p.λ), quantity
end
function _prepare_context_data(p::Union{SpectralBandPyrometer, TwoBandsRatioPyrometer}, quantity::TabularQuantity)
    l1, l2 = Tuple(p.λ)
    return subrange_view(l1, l2, quantity) # Возвращает кортеж или пару (_l, _e)
end
function _prepare_context_data(p::Union{SingleWavelengthPyrometer , TwoWavelengthRatioPyrometer}, 
                            quantity::TabularQuantity)
    _l = Tuple(p.λ)
    _e = get_single_wavelength_value(_l, nothing, quantity)
    return _l, _e
end
#SpectralQuantityPyrometricContext(p::P , quantity::TabularQuantity , i_measured::F) = TabularQuantityContext(p , quantity , i_measured)
const TabularQuantityContext{L , E , F , P} = SpectralQuantityPyrometricContext{L , E , F , P} where E <: AbstractVector
# SpectralBandPyrometer <=> ctx.e <: Planck.AbstractSpectralQuantity

function (ctx::SpectralQuantityPyrometricContext{L , E , F , P})(t) where {L <: NTuple{2 , D} , 
                                                            E <: AbstractSpectralQuantity , 
                                                            F , P <: SpectralBandPyrometer} where D <: Number 

    _to_halley(Planck.Dₜplanck_weighted(ctx.e , ctx.λ[1] ,ctx.λ[2], t) , ctx.i_measured) # 
end

function (ctx::GenericSpectralQuantityIntegrationContext{L , E , F , P})(t) where {L <: NTuple{2 , D} , 
                                                            E <: AbstractSpectralQuantity , 
                                                            F , P <: SpectralBandPyrometer} where D <: Number 

    #_to_halley(Planck.Dₜplanck_weighted(ctx.e , ctx.λ[1] ,ctx.λ[2], t) , ctx.i_measured) # 
end


# the same for all pyrometers tabular data
function (ctx::SpectralQuantityPyrometricContext{L , E})(t) where {L <: AbstractVector , 
                                                                    E <: AbstractVector}                                                            
    _to_halley(Planck.Dₜplanck_weighted(ctx.e , ctx.λ , t) , ctx.i_measured) # discrete integrator 
end

function (ctx::SpectralQuantityPyrometricContext{ <: Any , <: NTuple{1} , <:Any , P})(t) where { P <: SingleWavelengthPyrometer} 
    return _to_halley(
        Planck.Dₜibb(first(ctx.λ) , t) , 
        ctx.i_measured, 
        first(ctx.e)
    )
end
function (ctx::SpectralQuantityPyrometricContext{L, E , F  , P})(t) where {L , E <: AbstractSpectralQuantity , F , P <: SingleWavelengthPyrometer}
    l = first(ctx.λ)
    (i , di , ddi)  = Planck.Dₜibb(l , t)
    (e , de , dde) = Planck.eval_Dₜ(ctx.e , l , t)
    return _to_halley(
                    (
                        e * i, 
                        di * e + de * i , 
                        ddi * e + 2di*de + dde * i 
                    )
                    ,
                    ctx.i_measured)
end
# two-bands-ratio abstract continuous 
function (ctx::SpectralQuantityPyrometricContext{L , E})(t) where {L <: NTuple{2 , D}  , E <: AbstractSpectralQuantity} where D <: Tuple 
    tpl = Planck.Dₜplanck_weighted_ratio(ctx.e , ctx.λ[1] , ctx.λ[2] , t )
    return _to_halley(tpl, ctx.i_measured )
end
function (ctx::SpectralQuantityPyrometricContext{L , E})(t) where {L <: Tuple{D , D} , E <: Tuple{Q,Q} } where {D <: AbstractVector , Q<: AbstractVector}
    tpl = Planck.Dₜplanck_weighted_ratio(ctx.e[1] , ctx.λ[1] , ctx.e[2] , ctx.λ[2] , t )
    return _to_halley(tpl , ctx.i_measured )
end

function (ctx::SpectralQuantityPyrometricContext{L, E , F  , P})(t) where {L , E <: AbstractSpectralQuantity , F , P <: TwoWavelengthRatioPyrometer}
    l1 , l2 = ctx.λ[1],ctx.λ[2]
    
    (r , dr , ddr)  = Planck.Dₜspectral_ratio(l1 , l2 , t)
    (e1 , de1 , dde1) = Planck.eval_Dₜ(ctx.e , l1 , t)
    (e2 , de2 , dde2) = Planck.eval_Dₜ(ctx.e , l2 , t)
    
    slope = e1/e2
    dslope = Planck._spectral_ratio_first_derivative(e1 , de1 , e2 , de2)
    d2slope = Planck._spectral_ratio_second_derivative(e1 , de1 , dde1 , e2 , de2 , dde2)
    return _to_halley(
                        (
                            slope * r, 
                            dr * slope +  r * dslope, 
                            ddr * slope + 2 * dr * dslope + d2slope * r 
                        )
                    ,
                    ctx.i_measured)
end


# this version of wrapper is for the TwoBandsRatioPyrometer
function (ctx::SpectralQuantityPyrometricContext{L, E})(t) where {L <: Tuple{Number, Number}, E <: Tuple{Number, Number}}
    ratio_constant = ctx.e[1] / ctx.e[2]
    return _to_halley(Planck.Dₜspectral_ratio(ctx.λ[1], ctx.λ[2], t), ctx.i_measured, ratio_constant)
end

# Multiwavelength pyrometry 

measure(p::MultiWavelengthPyrometer , i::AbstractVector; 
                                starting_vector = nothing , kwargs...) = p.mwp(i; starting_vector = starting_vector)

#measure(p::MultiWavelengthPyrometer , i::Union{IsothermalSpectralQuantity , Abstract}; 
                        #starting_vector = nothing) = p.mwp(i.(wavelength(p)); starting_vector = starting_vector)



"""
    (p::AbstractPyrometer)(i; T_starting::Number=1000.0)

Various versions of pyrometers calling to get the measured temperature from external radiation 

Surface emissivity is taken from pyrometer 

`p(i::Number) ` - intensity provided as a single number 

`p(i::Union{IsothermalSpectralQuantity , AbstractDiscreteQuantity})` - intensity is provided as temperature independent discrete of continuos quantity 

Surface emissivity is provided externally as a continuous function [`AbstractContinuousOrDiscreteQuantity`](@ref)

`p(i::Number , ϵ::AbstractContinuousOrDiscreteQuantity)` - single number intensity 

`p(i::Union{AbstractDiscreteQuantity , IsothermalSpectralQuantity} , ϵ::AbstractContinuousOrDiscreteQuantity)` - continuous or discrete intensity (temperature imdependent)

Surface emissivity , measured intensity (temperature dependent with temperature `radiation_temperature`) is provided as a continuos or discrete quantity 

`p(i::AbstractContinuousOrDiscreteQuantity, 
    radiation_temperature::Number, 
    ϵ::AbstractContinuousOrDiscreteQuantity)` - continuous or discrete intensity (temperature imdependent)

    !Important 

    If the intensity is provided as a single number it should be the same quantity as the pyrometer acceps 
    e.g. if `p` is `TwoBandsRatioPyrometer` the intensity should be the ratio of two integral 
    (within pyrometer's working range) intensities. Example :
    ```julia
        import PlanckFunctions: band_power
        p = TwoBandsRatioPyrometer((2.0 , 3.0), (4.0 , 5.0))
        Ttrue = 1200.0
        i_number = band_power(1200 , λₗ = 2.0 , λᵣ = 3.0)/band_power(1200 , λₗ = 4.0 , λᵣ = 6.0)
        p(i_number) # returns Ttrue
        # Pyrometers.Planck (which is PlanckFunctions)
    ```

    On the opposite side , for continuous intensity (including temperature dependent) 
    it should be the intensity itself , the procedure of measurements 
    extracts pyrometer's specific signal from this function internally
      ```julia
        import PlanckFunctions: band_power , ibb
        p = TwoBandsRatioPyrometer((2.0 , 3.0) , (4.0 , 5.0))
        Ttrue = 1200.0
        i_iso = IsothermalSpectralQuantity(
                        Base.Fix2(ibb , Ttrue)
        ) # wrapper around planck spectral intensity (not spectral ratio)

        p(i_iso) # returns Ttrue
    ```  

"""
(p::AbstractPyrometer)(i; T_starting::Number=1000.0 ,  segbuf=nothing, kwargs...) = measure(p , i ;  T_starting = T_starting ,  segbuf=segbuf,  kwargs...)
(p::AbstractPyrometer)(i  , ϵ::Union{Number , NTuple{2}}; T_starting::Number=1000.0,  segbuf=nothing,  kwargs...) = measure(p , i  , ϵ ;  T_starting = T_starting,  segbuf=segbuf,  kwargs...)                    
(p::AbstractPyrometer)(imeasured, 
                        ϵ::AbstractContinuousOrDiscreteQuantity; T_starting = 600.0,  segbuf=nothing,  kwargs...) = measure(p , imeasured,  ϵ; T_starting = T_starting,  segbuf=segbuf,  kwargs...)

(p::AbstractPyrometer)(i::AbstractContinuousOrDiscreteQuantity , 
                    radiation_temperature::Number , 
                    ϵ::AbstractContinuousOrDiscreteQuantity; 
                    T_starting::Number = 1000.0 ,  segbuf=nothing,  kwargs...) = measure(p , i , radiation_temperature , ϵ ; T_starting = T_starting,  segbuf=segbuf,  kwargs...)