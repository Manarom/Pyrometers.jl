module PyrometersZygoteExt

    using Pyrometers, Zygote, ADTypes
    import Pyrometers.Planck: eval_Dₜ

    function eval_Dₜ(q::GenericDifferentiableSpectralQuantity{F, AutoZygote}, λ, T) where {F}
        val = q.f(λ, T)
        d1 = Zygote.gradient(t -> q.f(λ, t), T)[1]
        d2 = Zygote.hessian(t -> q.f(λ, t), T)[1] 
        return (val, d1, d2)
    end

end