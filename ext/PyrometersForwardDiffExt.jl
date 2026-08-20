module PyrometersForwardDiffExt

    using Pyrometers , ADTypes
    using ForwardDiff
    import Pyrometers: GenericDifferentiableSpectralQuantity
    import Pyrometers.Planck: eval_Dₜ

    
    function eval_Dₜ(q::GenericDifferentiableSpectralQuantity{Q , B}, λ, T) where {Q, B <: AutoForwardDiff}
        f_T(t) = q.f(λ, t)
        val = q.f(λ, T)
        d1  = ForwardDiff.derivative(f_T, T)
        d2  = ForwardDiff.derivative(t -> ForwardDiff.derivative(f_T, t), T)
        return (val, d1, d2)
    end

end # module