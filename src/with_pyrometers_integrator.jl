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
     intensity_function::AbstractSpectralQuantity; segbuf=nothing, kwargs...) 
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
function integrate(p::MultiWavelengthPyrometer{N} , intensity::Union{AbstractDiscreteQuantity , IsothermalSpectralQuantity}; kwargs...) where N 
    SVector{N}(intensity.(wavelengths(p)))
end 
function integrate(p::MultiWavelengthPyrometer{N} , t::Number , intensity::AbstractSpectralQuantity; kwargs...) where N 
    SVector{N}(intensity.(wavelengths(p) , t))
end    