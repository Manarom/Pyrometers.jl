module MultiWavelengthPyrometryExt

    using Optimization 
    import Pyrometers
    const OPTIMIZATION_FUN = Optimization.OptimizationFunction(Pyrometers.disc , 
                                                               grad =  Pyrometers.grad! , 
                                                               hess = Pyrometers.hess!)

    function Pyrometers._solve_problem(point ,
                                starting_vector , lb , ub , optimizer )

        probl= Optimization.OptimizationProblem(
                        OPTIMIZATION_FUN, 
                        starting_vector,
                        point,
                        lb=lb, 
                        ub=ub
                        )

        results = solve(probl , optimizer())
        return (results , problem)
    end    
end