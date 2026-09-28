module MultiWavelengthPyrometryOptimizationExt

    using Optimization , OptimizationOptimJL
    import Pyrometers
    const OPTIMIZATION_FUN = OptimizationFunction( Pyrometers.disc , 
                                                                grad =  Pyrometers.grad! , 
                                                                hess = Pyrometers.hess!
                                                               )
    const DEFAULT_OPTIMIZER = LBFGS
    Pyrometers._solve_problem(p::Union{Pyrometers.MWPPoint , Pyrometers.BBPoint} , sv , lb , ub , ::Nothing ) = Pyrometers._solve_problem(p , sv , lb , ub , DEFAULT_OPTIMIZER )
    
    function Pyrometers._solve_problem(point::Union{Pyrometers.MWPPoint , Pyrometers.BBPoint} ,
                                                starting_vector , 
                                                lb , ub , 
                                                optimizer )               

        probl = OptimizationProblem(
                            OPTIMIZATION_FUN, 
                            starting_vector,
                            point;
                            lb=lb, 
                            ub=ub
                        )

        results = solve(probl , optimizer())
        return (results , optimizer , probl)
    end    
end