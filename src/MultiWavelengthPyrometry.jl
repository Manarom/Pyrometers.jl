
module MultiwavelengthPyrometry
    using LinearAlgebra, #
    MKL, # using MKL turns default LinearAlgebra from library from openBLAS to mkl  
    Optimization,
    OptimizationOptimJL, 
    #Interpolations,
    StaticArrays,
    RecipesBase,
    Distributions,
    .Pyrometers
    
    #Base.convert(::Type{SciMLBase.ReturnCode.T}, s::Symbol) = SciMLBase.ReturnCode.T(s)
    using ScaledPolynomials

    import PlanckFunctions as Planck
    include("MultiwavelengthPyrometryTypes.jl") # Brings types and functions for working with types
    #include("Pyrometers.jl") 
    const NPOINT = 30

    export  MWPPoint,# type for least-square fitting 
            EmPoint, # type for BB temperature fitting
            fit_T!,# function to fit the BB and real surface temperature
            Pyrometers, # pyrometers module
            optimizer_switch, # supported optimizer selection function (returns optimizer constructor)
            optim_dic,
            support_constraint_optimizers,
            support_lagrange_constraints,
            em_cons!, # emissivity constraint function (returns maximum and minimum values of the emissivity in the whole spectral range)
            emissivity!, # fills the emissivity
            feval!,
            jacobian!,
            hess!,
            hess_approx!,
            grad!,
            residual!,
            disc
    """
    All supported optimizers
    """
    const optim_dic = Base.ImmutableDict("NelderMead"=>NelderMead,
                            "Newton"=>Newton,
                            "BFGS"=>BFGS,
                            "GradientDescent"=>GradientDescent,
                            "NewtonTrustRegion"=>NewtonTrustRegion,
                            "ParticleSwarm"=>ParticleSwarm,
                            "Default"=>NelderMead,
                            "LBFGS"=>LBFGS,
                            "IPNewton"=>IPNewton) # list of supported optimizers
    
    const DEFAULT_OPTIMIZER = Ref(LBFGS())
    
    """
    Optimizers supporting box-constraint optimization
    """                        
    const support_constraint_optimizers = ["NelderMead", 
                                            "LBFGS",
                                            "IPNewton",
                                            "ParticleSwarm"]
    """
    Optimizer supporting lagrangian-constraint optimization
    """                                        
    const support_lagrange_constraints = ["IPNewton"]
    sumabs2(a) = sum(abs2,a)
    """
    optimizer_switch(name::String;is_constraint::Bool=false,
                is_lagrange_constraint::Bool=false)

Returns the appropriate optimizer constructor
Input:
    name - the name of th eoptimizer, 
    is_constraint - is box constraint problem formulation, 
    is_lagrange_constraint - use Lagrange constraints (supported only by IPNewton)

Returns the appropriate optimizer constructor
Input:
    name - the name of th eoptimizer, 
    is_constraint - is box constraint problem formulation, 
    is_lagrange_constraint - use Lagrange constraints (supported only by IPNewton)
"""
    function optimizer_switch(name::String;is_box_constraint::Bool=false,
                is_lagrange_constraint::Bool=false)
            if is_lagrange_constraint
                name = filter(x -> ==(name,x),support_lagrange_constraints)
                default_ = optim_dic[support_lagrange_constraints[1]]
                return length(name)>=1 ? get(optim_dic,name[1],default_ ) :
                                                                    default_            
            elseif is_box_constraint
                name = filter(x -> ==(name,x),support_constraint_optimizers)
                return length(name)>=1 ? get(optim_dic,name[1],optim_dic["Default"]) :
                optim_dic["Default"]
            else
                return get(optim_dic,name,optim_dic["Default"])
            end
    end
    ## BAND PYROMETRY POINT METHODS


    const LAGRANGE_OPTIM_FUN = OptimizationFunction(disc,grad=grad!,hess=hess!,cons=em_cons!)
    const OPTIM_FUN  = OptimizationFunction(disc,grad=grad!,hess=hess!) 
    """
    fit_T!(point::Union{EmPoint,MWPPoint};
            optimizer_name::String="Default",
            is_box_constraint::Bool=false,
            is_lagrange_constraint::Bool=false, 
            emissivity_range::C=nothing, 
            temperature_range::B=nothing) where {B <: Union{AbstractVector,Nothing,NTuple{2}},C <: Union{AbstractVector,Nothing,NTuple{2}}}

    Input:
        point - (modified) real suraface (MWPPoint) of blackbody thermal emission object
        (optional)
        optimizer_name - the name of optimizer must be the key of optim_dic
        is_constraint - is box-constraint flag
        is_lagrange_constraint - is Lagrange-constraint flag
    Returns:
        if the point type is EmPoint, the output is:
            named tuple with (T - fitted temperature,
                                res - optimization output object,
                                optimizer - chosen optimizer)
        if the point is of the MWPPoint type, the output is:
            named tuple with (T - fitted temperature ,a - fitted emissivity approximation coefficients,
                                ϵ - emissivity spectrum in the whole wavelemgth range,
                                res - optimization output object,
                                optimizer - chosen optimizer)                
"""
function fit_T!(point::Union{EmPoint , MWPPoint};
            optimizer_name::String="Default",
            is_box_constraint::Bool=false,
            is_lagrange_constraint::Bool=false, 
            emissivity_range::C=nothing, 
            temperature_range::B=nothing , 
            kwargs...) where {B <: Union{AbstractVector,Nothing,NTuple{2}},C <: Union{AbstractVector,Nothing,NTuple{2}}}
            
            !(optimizer_name == "Default")  || return fitT_default(point)


            optimizer = optimizer_switch(optimizer_name,
                                    is_box_constraint = is_box_constraint,
                                    is_lagrange_constraint = is_lagrange_constraint)
        fun = is_lagrange_constraint ? LAGRANGE_OPTIM_FUN : OPTIM_FUN
        # by default all derivatives are supported
        if point isa EmPoint
            starting_vector = MVector{1}([235.0])
        else
            starting_vector = copy(point.x);
        end
        if is_lagrange_constraint
            (lb , ub) = evaluate_lagrange_constraints(point,emissivity_range,temperature_range)
            probl= OptimizationProblem(fun, 
                            starting_vector,
                            point, 
                            lcons = lb, # both min and max of emissivity should be not smaller than zero
                            ucons = ub, kwargs...) # both min and max should be higher than one        
        elseif is_box_constraint
            (lb , ub) = evaluate_box_constraints(point, emissivity_range, temperature_range)
            trim_starting_vector_to_box!(starting_vector , lb , ub)
            probl= OptimizationProblem(fun, 
                            starting_vector,
                            point,
                            lb=lb, 
                            ub=ub, kwargs...)
        else # unconstraint
            probl= OptimizationProblem(fun, 
                                starting_vector,
                                point; kwargs...)           
        end
        results = solve(probl,optimizer())
        isa(point, MWPPoint) ? copyto!(point.x , results.u) : feval!(point,results.u)
        return  fitting_result(point, results, optimizer) 
                        
    end
    
   
    function fitT_default(point::EmPoint)
        probl= OptimizationProblem(OPTIM_FUN, MVector{1}([235.0]),point)
        results = solve(probl,DEFAULT_OPTIMIZER[])   
        feval!(point,results.u)
        return (T=temperature(point), res=results, optimizer=DEFAULT_OPTIMIZER[])                 
    end
    function fitT_default(point::MWPPoint)
        probl= OptimizationProblem(OPTIM_FUN, point.x,point)
        results = solve(probl,DEFAULT_OPTIMIZER[])   
        feval!(point,results.u)
        return (T=temperature(point), res=results, optimizer=DEFAULT_OPTIMIZER[])                 
    end
    function (emp::EmPoint)(I::AbstractVector) 
        copyto!(emp.I_measured , I)
        return fitT_default(emp).T
    end
    function (emp::EmPoint)(eps::Number) 
        emp.I_measured ./= eps # dividing by emissivity 
        # copyto!(emp.I_measured,I)
        T =  fitT_default(emp).T
        emp.I_measured .*= eps
        return T
    end
    function (emp::MWPPoint)(I::AbstractVector)
        copyto!(emp.e_p.I_measured,I)
        return fitT_default(emp).T
    end
    function (emp::Union{EmPoint,MWPPoint})()
        return fit_T!(emp).T
    end
    """
    covariance(bp::MWPPoint)

Evaluates the covariance matrix as Cov(x) = 2σ²H⁻¹
"""
function fitting_covariance(bp::MWPPoint{N,Nx3,P}) where {N,Nx3,P}
    sigma_square = sumabs2(bp.e_p.ri)/degrees_of_freedom(bp)
    h = similar(bp.hessian)
    hess!(h, bp.x, bp::MWPPoint)
    return 2*sigma_square*inv(h)
end
"""
    fitting_covariance(em::EmPoint{N})

Evaluates the covariance matrix as Cov(x) = 2σ²H⁻¹
"""
function fitting_covariance(em::EmPoint{N,Nx3,T}) where {N,Nx3,T}
    sigma_square = sumabs2(em.ri)/degrees_of_freedom(em)
    h = MMatrix{1,1,T,1}(undef)
    hess!(h,temperature(em),em)
    return 2*sigma_square*inv.(h)
end
fitting_variance(em::EmPoint) = vec(fitting_covariance(em))
"""
    fitting_variance(bp::MWPPoint)

Returns the optimization variable variance (diagonal of the covariance matrix)
"""
fitting_variance(bp::MWPPoint) = collect(diag(fitting_covariance(bp)))
fitting_error(p::Union{MWPPoint,EmPoint};probability = 0.95,only_std::Bool=false) =(only_std ? 1.0 : student_coefficient(degrees_of_freedom(p),probability))*sqrt.(fitting_variance(p))
"""
    student_coefficient(degrees_of_freedom::Int, probability; digits::Int = 3, side::Int = 2)

Evaluates Student's distribution coefficient
"""
function student_coefficient(degrees_of_freedom::Int, probability;  side::Int = 2)
	if side == 2
        probability = (1 + probability)/2
    end
	return Distributions.quantile(Distributions.TDist(degrees_of_freedom), probability)
end
    @recipe function f(m::EmPoint)
        minorgrid--> true
        gridlinewidth-->2
        dpi-->600
        xlabel-->"Wavelength"
        ylabel-->"Spectral intensity"
        linewidth-->3

        @series begin 
            label:= "Measured"    
            linewidth:=2
            markershape:=:none
            fillrange:=0
            fillalpha:=0.3
            (m.λ,m.I_measured)
        end
        @series begin 
            label:="Fitted"    
            linewidth:=2
            fillrange:=0
            fillalpha:=0.3
            markersize := 3
            markershape:=:diamond
            (m.λ,m.Ib)
        end
    end
    @recipe function f(m::MWPPoint)
        minorgrid--> true
        gridlinewidth-->2
        dpi-->600
        xlabel-->"Wavelength"
        ylabel-->"Spectral intensity"
        linewidth-->3

        @series begin 
            label:= "Measured"    
            linewidth:=2
            markershape:=:none
            fillrange:=0
            fillalpha:=0.3
            (m.e_p.λ,m.e_p.I_measured)
        end
        @series begin 
            label:="Fitted"    
            linewidth:=2
            fillrange:=0
            fillalpha:=0.3
            markersize := 3
            markershape:=:diamond
            (m.e_p.λ,m.Ic)
        end
    end
end

