module PyrometersDataInterpolationsExt

    using DataInterpolations, Pyrometers
    import Pyrometers.QuadGK

    function Pyrometers.QuadGK.quadgk(f::IsothermalSpectralQuantity{SQ} , 
        a::Number, b::Number ;kwargs...) where SQ <: DataInterpolations.AbstractInterpolation{T} where T

        X = eltype( f.f.t )    
        val = DataInterpolations.integral(f.f , convert(X , a), convert(X , b))
        return (val, zero(val))       
    end

end
