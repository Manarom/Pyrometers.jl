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

