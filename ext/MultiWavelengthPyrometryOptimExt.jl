module MultiWavelengthPyrometryOptimExt

    using Optim
    using StaticArrays
    using Pyrometers
    using NLSolversBase
    import Pyrometers as P

    const DEFAULT_MWOPTIMIZER = LBFGS
    const DEFAULT_BBOPTIMIZER = Brent
    function P._solve_problem(p::P.BBPoint , _ , lb::Number , ub::Number , ::Nothing ) 
            f_scalar(x) = Pyrometers.disc(SVector{1}(x), p)
            return (optimize(f_scalar, lb , ub , Brent()) , Brent , nothing)
    end
    P._solve_problem(p::P.BBPoint , _ , 
                        lb::AbstractVector , 
                        ub::AbstractVector , 
                        ::Nothing )  = P._solve_problem(p , nothing , last(lb) , last(ub) , nothing )
    function P._solve_problem(point::Union{P.BBPoint , P.MWPPoint} ,
                                                starting_vector , 
                                                lb , 
                                                ub , 
                                                optimizer )               

            f = NLSolversBase.only_fgh!(make_fgh!(point))
            if lb !== nothing && ub !== nothing
                result = optimize(f ,   lb , ub , starting_vector, Fminbox(optimizer()))
            else
                result = optimize(f,  starting_vector , optimizer())
            end
        return (result , optimizer , nothing)
    end 
    function make_fgh!(point)
        return function fgh!(F , G , H , x)

            isnothing(G) || Pyrometers.grad!(G , x , point)
            isnothing(H) || Pyrometers.hess!(H , x , point)
            !isnothing(F) && return Pyrometers.disc(x, point)
            return nothing
        end

    end
    Pyrometers._solve_problem(p::P.MWPPoint , sv , lb , ub , ::Nothing ) = Pyrometers._solve_problem(p , sv , lb , ub , DEFAULT_MWOPTIMIZER )
end