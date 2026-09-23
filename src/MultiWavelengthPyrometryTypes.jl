
# MultiwavlengthPyrometryTypes should be included in the MultiwavlengthPyrometry module
const DEFAULT_TEMPERATURE_RANGE = Ref((20.0 , 10000.0)) # default temperture range used for bounded optimization
const DEFAULT_EMISSIVITY_RANGE = Ref((0.1 , 1.0))
struct DefaultOptimizer end
(d::DefaultOptimizer)() = d
#function 
"""
BBPoint type stores data on thermal emission spectrum and its 
first and second derivatives it also stores "Measurements " vector 
which further can be fitted? it also provides the constructor
BBPoint(I_measured,λ) -  I_measured is a measured spectrum
                      -  λ - wavelength vector (in μm)

"""
struct BBPoint{N , Nx3 , T} 
    I_measured::MVector{N,T}# data to fit
    λ:: MVector{N,T}  # wavelength vector (it is constant during the optimization)
    Ib::MVector{N,T} # Planck function values vector ????
    ri::MVector{N,T} # residual vector
    r::Base.RefValue{T} # discrepancy value
    ∇I::MVector{N,T} # first derivative value
    ∇²I::MVector{N,T} # second derivative vector
    amat::MMatrix{N,3,T,Nx3} # intermediate private data used to speed up the planck function evaluation
    # temperatures of:
    Tib::Base.RefValue{T} # Planck intensity  
    Tri::Base.RefValue{T} # Residual vector  
    T∇ib::Base.RefValue{T} # Planck derivative  
    Tgrad::Base.RefValue{T} # Gradient of emission discrepancy function 
    T∇²ib::Base.RefValue{T} # Planck function second derivative 
    Thess::Base.RefValue{T} # Discrepancy function second derivative
    """
    BBPoint(I_measured::AbstractVector,λ::AbstractVector)

Constructor of the BBPoint object instance
Input: 
    I_measured - mesured blackbody spectral intensity
    λ - wavelength in μm
"""
function BBPoint(I_measured::StaticArray{Tuple{N},T,1},λ::StaticArray{Tuple{N},T,1}) where {N,T}

       #points_number = length(λ)
       SVectType = SVector{N,T}
       MVectType = MVector{N,T}
       MatType = MMatrix{N,3,T,3*N} 

       new{N,3*N,T}(
            SVectType(I_measured), # measured value
            MVectType(λ),# wavelength
            MVectType(undef),#Ib::AbstractVector# Planck function values vector
            MVectType(undef),#ri::AbstractVector  # discrepancy vector
            Ref(maxintfloat(Float64)),#r::Base.RefValue{Float64}# discrepancy value
            MVectType(undef),#∇I::AbstractVector # first derivative value
            MVectType(undef),#∇²I::AbstractVector # second derivative vector
            MatType(zeros(N,3)),#amat::Matrix{Float64} # intermediate private data 
            Ref(0.0),#  stores the temperature of Planck function evaluation
            Ref(0.0), # Tri
            Ref(0.0),# T∇ibb
            Ref(0.0),# Tgrad
            Ref(0.0),# T∇²ibb
            Ref(0.0),# Tsec          
       ) # calling the constructor

    end
end
BBPoint(i::AbstractVector , l::AbstractVector) = begin 
    N = length(i)
    (length(i) == length(l)) || error("vectors must be of the same length") 
    return BBPoint(SVector{N}(i) , SVector{N}(l))
end
#VanderMatrix(em::BBPoint,vv::Val{CN};poly_type::Symbol = :stand) where CN = VanderMatrix(em.λ , vv ,poly_type = poly_type)

"""
    MWPPoint type stores data of thermal emission spectrum of a real body with 
emissivity polynomial approximation, and  its first and second derivatives
it also stores "measurements" vector which further can be fitted, it also 
provides the constructor MWPPoint(I_measured,λ,initial_x,polynomial_type) 
where:
    -  I_measured is a measured spectrum
    -  λ - wavelength vector (in μm)
    - initial_x - starting parameters vector (initial_x[end] - starting temperature,
        x[1:end-1] - emissivity approximation coefficients)
    - polynomial_type - string of polynomial (this value governs the Vandermonde matrix form)

N - number ob wavelength points
Nx3 - 3*N used to store the intermediate data for planck function and derivarives evaluation
L - length of optimization variables vector
P -   number of variables to be optimized (can be less or equal to the number of )
NxP - N*P number of elements in jacobian 
PxP - P*P number of elements in hessian
Pm1 - P-1 number of parameters approximating 
NxPm1 - N*(P-1) number of vandermatrix elements
T - type of data
"""
struct MWPPoint{N , Nx3 ,P, NxP, PxP, Pm1 , NxPm1, Pm1xPm1 , T , PolyType}#{N,P,T} # N - wavelength number, CN - parameters number + 1
    # N , Nx3 , P, NxP, PxP, Pm1 , NxCN, CNxCN , T
    # Stores data about the spectral band, BBemission spectrum and experimental measured spectrum
    bb::BBPoint{N,Nx3,T} 
    # Additional data storages
    x::MVector{P,T}   #Px1 # Optimization variables vector
    # x[end] - temperature, x[1:end-1] - emissivity poynomial approximation
    Ic::MVector{N,T} #Lx1 # sample spectral emittance
    Iₛᵤᵣ::SVector{N,T} #::Lx1 # surrounding radiation spectra
    ϵ::MVector{N,T}  #Lx1 # spectral emissivity in band
    r::MVector{N,T} #Lx1 # residual vector
    jacobian::MMatrix{N,P,T,NxP}#LxP # Jacobian matrix
    hessian_approx::MMatrix{P,P,T,PxP} #PxP # approximate hessian matrix
    hessian::MMatrix{P,P,T,PxP} # Hesse matrix
    vandermonde::VanderMatrix{N, Pm1, T, NxPm1, Pm1xPm1 , PolyType} # Vandermonde matrix type VanderMatrix{N,CN,T,NxCN,CNxCN} - N - rows number, CN - columns of vander number ()
    # internal usage
    x_em_vec::MVector{P,T} # vector of emissivity evaluation values
    x_jac_vec::MVector{P,T} #Px1 # vector of jacobian claculation parameters
    x_hess_approx::MVector{P,T}#Px1 # vector of the approximate hessian calculation
    x_hess_vec::MVector{P,T} #Px1 # vector of hessian calculation parameters
    is_has_Iₛᵤᵣ::Bool # flag 

"""
    MWPPoint(measured_Intensity::AbstractVector,
                        λ::AbstractVector,
                        initial_x::AbstractVector;
                        polynomial_type::Symbol="stand",
                        I_sur::AbstractVector=[])

Constructor for band pyrometry fitting, 
    λ - wavelength vector, 
    initial_x - starting optimization vector 
    polynomial_type - type of polynomial for emissivity approximation
"""
function MWPPoint(measured_Intensity::StaticArray{Tuple{N},T,1},
                        λ::StaticArray{Tuple{N},T,1},
                        initial_emissivity::StaticArray{Tuple{Pm1} , T , 1},
                        initial_temperature::T,
                        ::Type{PolyType} = BernsteinSymPoly{Pm1,T};
                        I_sur::Union{StaticArray{Tuple{N},T,1},Nothing}=nothing) where PolyType <: AbstractPoly{Pm1,T} where {N,Pm1,T}

       #PolyTypeAbs = haskey(ScaledPolynomials.SUPPORTED_POLYNOMIAL_TYPES , polynomial_type) ? ScaledPolynomials.SUPPORTED_POLYNOMIAL_TYPES[polynomial_type] : BernsteinSymPoly

       # if entered polynomial type is not supported then it turns to "simple"
       #L = length(λ) #total number of spectral points
       #P = length(initial_x) # full number of the optimization variables
       #polynomial_degree =  P - 2 #degree of emissivity polynomial approximation
       # polynomial degree goes from 0,1... where 1 is linear approximation
       # {N , Nx3 , P, NxP, PxP, Pm1 , NxPm1, Pm1xPm1 , T}
       initial_x = (initial_emissivity... , initial_temperature)
       P = Pm1 + 1
       Nx3 = 3*N
       NxP = N*P
       PxP = P*P 
       #Pm1 = P - 1
       NxPm1 = N*(P - 1)
       Pm1xPm1 = (P - 1)*(P - 1) 
       
       # @show isconcretetype(PolyTypeAbs{Pm1})

       Nx1_T = MVector{N,T} # independent data column
       Px1_T = MVector{P,T} # optimization variables vector 
       NxP_T = MMatrix{N,P,T,NxP} # Jacobian type
       PxP_T = MMatrix{P,P,T,PxP} # Hessian type     
       #LxPm1_T = MMatrix{N,Pm1,T,NxPm1} #Vandermonde matrix type
       is_has_Iₛᵤᵣ = !isnothing(I_sur) && length(I_sur)==length(λ)
       Isr =  is_has_Iₛᵤᵣ ? SVector{N}(I_sur) : SVector{N}(zeros(T,N))
       # {N , Nx3 , P, NxP, PxP, Pm1 , NxPm1, Pm1xPm1 , T}
       #PolyTypeAbs = get(ScaledPolynomials.SUPPORTED_POLYNOMIAL_TYPES, polynomial_type, ScaledPolynomials.BernsteinSymPoly)
       poly = PolyType()
       #PolyType = PolyTypeAbs{Pm1 ,T}
       new{N , Nx3 , P, NxP, PxP, Pm1 , NxPm1, Pm1xPm1 , T , PolyType}(
                BBPoint(measured_Intensity,λ),# filling BB emission obj
                Px1_T(initial_x), #em_poly
                Nx1_T(undef), # emissivity
                Nx1_T(undef), # Ic corrected emission spectrum
                Isr, # Iₛᵤᵣ surrounding radiation exclusion
                Nx1_T(undef), # r,residual vector function
                NxP_T(undef),# jacobian
                PxP_T(undef),# approximate hessian
                PxP_T(undef),# hessian
                VanderMatrix(SVector{N}(λ), # wavelength
                            poly 
                ),
                Px1_T(undef), # x_em_vec - emissivity evaluation vector
                Px1_T(undef), # x_jac_vec - Jacobian evaluation vector
                Px1_T(undef), # x_hess_approx - approximate Hessian evaluation vector
                Px1_T(undef), # x_hess_vec - rigorous hessian evaluation vector
                is_has_Iₛᵤᵣ  # is_has_Iₛᵤᵣ - flag is true if the point has surrounding radiation correction part
                )
    end
end

temperature(emp::BBPoint) = emp.Tib[]
temperature(bp::MWPPoint) = bp.bb.Tib[]

pointsnumber(::Union{BBPoint{N},MWPPoint{N}}) where N = N
parnumber(::MWPPoint{N, Nx3, P}) where  {N, Nx3, P} = P
parnumber(::BBPoint) = 1
degrees_of_freedom(p::Union{BBPoint,MWPPoint}) = pointsnumber(p) - parnumber(p)
emissivity(p::MWPPoint) = copy(p.ϵ)

#function emissivity_polynomial(p::MWPPoint{})

function clear_cache!(p::BBPoint{N , M , T}) where {N,M,T}

    protected = (:I_measured, :λ) 

    for field in propertynames(p)
        field in protected && continue
        
        val = getproperty(p, field)
        if val isa StaticArray
            val .= zero(T)          # In-place reset for MVectors/MMatrices
        elseif val isa Base.RefValue
            val[] = zero(T)         # In-place reset for Refs
        end
    end
    residual!(p , first(DEFAULT_TEMPERATURE_RANGE[]))
    return p
    
end
const SimplyTypedMW{N,P,T} = MWPPoint{N, Nx3, P, NxP, PxP, Pm1, NxPm1, Pm1xPm1, T} where {N, Nx3, P, NxP, PxP, Pm1, NxPm1, Pm1xPm1, T}
function clear_cache!(p::SimplyTypedMW{N , P , T}) where {N , P , T}
    
    clear_cache!(p.bb)

    protected = (:Iₛᵤᵣ, :is_has_Iₛᵤᵣ ,:vandermonde, :bb) 

    for field in propertynames(p)
        field in protected && continue
        
        val = getproperty(p, field)
        val .= zero(T) 
    end
    @. p.r = p.bb.ri
    return p
    
end
#emissivity(p::BandPyrometryPoint,λ::AbstractVector) = 
    """
    box_constraints(bp::MWPPoint)

Evaluates box-constraint of the problem
"""
function evaluate_box_constraints(bp::MWPPoint{N, Nx3, P, NxP, PxP, Pm1, NxPm1, Pm1xPm1, T},
        emissivity_range::B  = nothing,
        temperature_range::C = nothing )  where {B <: Union{Nothing , NTuple{2 , T} ,  NTuple{2 , <: Union{NTuple{Pm1 ,T} , StaticVector{Pm1 , T}} } } ,
                                                C <: Union{Nothing , NTuple{2,T}} } where {N, Nx3, P, NxP, 
                                                            PxP, Pm1, NxPm1, Pm1xPm1, T}
    # method calculates box constraints 
    # of the feasible region (dumb version)
        lb = copy(bp.x)
        ub = copy(bp.x)
        e_lb = @view lb[1:end - 1] 
        e_ub = @view ub[1:end - 1] 
        # b_all = isnothing(emissivity_range) ? (0.0 , 1.0) : (first(emissivity_range), last(emissivity_range))
        fill_emissivity_box_constraint!(e_lb , e_ub , bp.vandermonde, emissivity_range)
        # @show e_lb , e_ub 
        (lb[end], ub[end]) = if  isnothing(temperature_range) 
            extract_temperature_range(bp , emissivity_range) 
        else 
            (first(temperature_range) , last(temperature_range))
        end
    return (lb=lb , ub=ub)
end
extract_temperature_range(::MWPPoint , ::Nothing) = DEFAULT_TEMPERATURE_RANGE[]
function extract_temperature_range(p::MWPPoint , e_range::NTuple{2 , <: Union{NTuple , StaticVector} } )
    lb = first(e_range)
    ub = last(e_range)
    return extract_temperature_range(p , (minimum(lb) , maximum(ub)))
end
function extract_temperature_range(p::MWPPoint , emissivity_range::NTuple{2,T}) where T
    (e1 , e2) = emissivity_range
    (e1 > e2) && ((e1 , e2) = (e2 , e1)) 
    return ( try_emissivity(p.bb , e2) , try_emissivity(p.bb , e1) ) # the lower and the upper limits on temperature
end

  """
    em_cons!(constraint_value::AbstractArray,
                            x::AbstractVector, 
                            bp::MWPPoint)


In-place filling of two-elemnt vector of [minimum,maximum] emissivity in the whole 
wavelength range  
This function is used in the constraints
Inputs:
    constraint_value - (modified)  two-element vector to be modified in-place
    x - optimization variables vector, x=[a1...an,T], where a1...an - emissivity approximations coefficients, T  - temperature 
    bp - (modified) 
"""
function em_cons!(constraint_value::AbstractArray,
                        x::AbstractVector, 
                        bp::MWPPoint)
    # evaluate the constraints on emissivity (it should not be greater than one in a whole spectra range)
    feval!(bp,x)  
    constraint_value .= extrema(bp.ϵ) # (minimum,maximum) values of the emissivity 
    return constraint_value
    #   in a whole spectrum range
end
    """
    emissivity!(bp::MWPPoint,x::AbstractVector)


Fills emissivity for the current BandPyrometry point
Input:
    bp - (modified) current spectral band pytometry point
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature   
Fills emissivity for the current BandPyrometry point
Input:
    bp - (modified) current spectral band pytometry point
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature   
"""
function emissivity!(bp::MWPPoint , x::AbstractVector)
    a = @view x[1:end-1] #emissivity approximation variables
    return mul!(bp.ϵ , bp.vandermonde.v , a)
end


    
    """
    feval!(bp::MWPPoint,x::AbstractVector)

Fills both the emissivity and the thermal emission spectrum for the current BandPyrometry point

Input:
    bp - (modified) current spectral band pytometry point
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature
"""
    function feval!(bp::MWPPoint , x::AbstractVector)
        # evaluates residual vector
        #a = @view x[1:end-1] #emissivity approximation variables
        feval!(bp.bb , x[end]) # refreshes planck function values
        if x != bp.x_em_vec # x_em_vec - emissivity calculation vector
            emissivity!(bp,x)
            if bp.is_has_Iₛᵤᵣ # has surrounding radiation correction
                @. bp.Ic = (bp.bb.Ib - bp.bb.Iₛᵤᵣ) * bp.ϵ # I=(Ibb-Isur)*ϵ
            else
                @. bp.Ic = bp.bb.Ib * bp.ϵ # I=Ibb*ϵ
            end
            bp.x_em_vec .= x
        end
        return bp.Ic
    end
        """
        residual!(bp::MWPPoint,x::AbstractVector)

    Fills emissivity, thermal emission spectrum and evaluates the residuals vector 
    for the current BandPyrometry point

    Input:
        bp - (modified) current spectral band pytometry point
        x - optimization variables vector, x=[a1...an,T],
        where a1...an - emissivity approximations coefficients, T  - temperature    
    """
    function residual!(bp::MWPPoint , x::AbstractVector)
        feval!(bp , x)   # feval! calculates function value only if current x is not the same as 
        @. bp.r =bp.bb.I_measured - bp.Ic # measured data - calculated 
        bp.bb.r[] = 0.5 * norm(bp.r)^2 # discrepancy value
        return bp.r # returns residual vector
    end
    #function residual!(p::BBPoint)
    """
    disc(x::AbstractVector,bp::MWPPoint)

Fills emissivity, thermal emission spectrum,evaluates the residuals vector
and calculates its norm for the current BandPyrometry point
Input:
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature 
    bp - (modified) current spectral band pytometry point
Fills emissivity, thermal emission spectrum,evaluates the residuals vector
and calculates its norm for the current BandPyrometry point
Input:
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature 
    bp - (modified) current spectral band pytometry point
"""
    function  disc(x::AbstractVector,bp::MWPPoint)
        residual!(bp,x)
        return bp.bb.r[]# returns current value of discrepancy
    end

"""
    jacobian!(bp::MWPPoint , x::AbstractVector)


Fills the Jacobian matrix for current bandpyrometry point
Input:
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature
    bp - (modified, stores Jacobian internally) current spectral band pytometry point 
Fills the Jacobian matrix for current bandpyrometry point
Input:
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature
    bp - (modified, stores Jacobian internally) current spectral band pytometry point 
"""
function jacobian!(bp::MWPPoint , x::AbstractVector) # evaluates Planck function
        ∇!(bp.bb , x[end]) # refresh Planck function first derivative
        if x != bp.x_jac_vec
            J1 = @view bp.jacobian[: , 1:end-1] # Jacobian without temperature derivatives
            J2 = @view bp.jacobian[:,end] # Last column of the jacobian 
            #a  = @view (x,1,end-1)
            eps_vec = emissivity!(bp, x)
            @. J1 = bp.bb.Ib * bp.vandermonde.v # diag(ibb)*V
            @. J2 = bp.bb.∇I * eps_vec# 
            copyto!(bp.x_jac_vec , x) # refresh jacobian calculation vector
        end
        return bp.jacobian
    end   

    """
    grad!(g::AbstractVector , x::AbstractVector , bp::MWPPoint)

In-place filling of the gradient vector of MWPPoint at point x
Input:
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature
    bp - (modified, recalculates residual vector and Jacobian if the 
    currently stored value was obtaibed for another optimization variables array)
    current spectral band pytometry point     
In-place filling of the gradient vector of MWPPoint at point x
Input:
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature
    bp - (modified, recalculates residual vector and Jacobian if the 
    currently stored value was obtaibed for another optimization variables array)
    current spectral band pytometry point     
"""
    function grad!(g::AbstractVector , x::AbstractVector  , bp::MWPPoint)
        residual!( bp , x)
        jacobian!( bp , x) # calculated Jₘ
        g .= - transpose(bp.jacobian) * bp.r # calculates gradient ∇f = -Jₘᵀ*r
        # mul!(g, transpose(bp.jacobian), bp.r, -1.0, 0.0) # gives an error ! (bug ?)
        #J_static = SMatrix(bp.jacobian)
        #r_static = SVector(bp.r)
        
        #g .= - (transpose(J_static) * r_static)
        return g
    end

    """
    hess_approx!(ha, x::AbstractVector,bp::MWPPoint)

In-place filling of the approximate hessian (Hₐ = Jᵀ*J (J - Jacobian)) 
of MWPPoint at point vector x, approximate Hessian can be used 
in optimization methods to approximate the full Hessian (e.g. in Gauss-Newton
or Levenberg-Marquardt methods)

Input:
    ha - Hessian matrix to be filled in-place
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature
    bp - (modified) current spectral band pytometry point     
"""
function hess_approx!(ha , x::AbstractVector  , bp::MWPPoint)
        # calculates approximate hessian which is Hₐ = Jᵀ*J (J - Jacobian)
        ha .= hess_approx!(bp, x) 
        return ha
    end
    function hess_approx!(bp::MWPPoint , x::AbstractVector )
        # calculates approximate hessian which is Hₐ = Jᵀ*J (J - Jacobian)
        if x != bp.x_hess_approx
            jacobian!(bp , x)
            bp.hessian_approx .= transpose(bp.jacobian) * bp.jacobian 
            # this matrix is always symmetric positive definite
            copyto!(bp.x_hess_approx , x)
        end
        return bp.hessian_approx
    end    
    """
    hess!(h , x::AbstractVector,bp::MWPPoint)
    
In-place filling of the whole hessian matrix for MWPPoint at point x
Input:
    ha - Hessian matrix to be filled in-place
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature
    bp - (modified) current spectral band pytometry point      
In-place filling of the whole hessian matrix for MWPPoint at point x
Input:
    ha - Hessian matrix to be filled in-place
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature
    bp - (modified) current spectral band pytometry point      
"""
function hess!(h , x::AbstractVector , bp::MWPPoint)
        hess!(bp , x)
        copyto!(h , bp.hessian) # filling external matrix with internally stored hessian
        return h
    end

    function hess!(bp::MWPPoint , x)
        if x != bp.x_hess_vec
            hess_approx!(bp.hessian , x , bp) # refresh the approximate hessian 
            # and fill hessian with approximate hessian Jᵀ*J
            # refreshes second derivative of the Planck function
            ∇²!(bp.bb , x[end]) 
            # H = Ha - Hm, Ha is approximate Hessian
            # Hm_vec = Vᵀ*I'ᴰ*r - vector Hm,
            # V - Vandermonde matrix, I'ᴰ - first 
            # derivative diagonal matrix,
            # r - residual vector
            last_hess_col = @view bp.hessian[1:end-1, end] 
            # view of the last column of the hessian 
            # initial formula: Hm_vec = Vᵀ*I'ᴰ*r  => transpose(V)*diagm(I')*r 
            # A*diagm(b) <=> A.*transpose(b) <=> transpose(Aᵀ.*b) 
            # Hm_vec = (V.*I')ᵀ*r
            last_hess_col .-= transpose(bp.vandermonde.v .* bp.bb.∇I) * bp.r
            bp.hessian[end , 1:end-1] .= last_hess_col # the sample
            # only right-down corner of hessian contains the second derivative
            # hm = rᵀ*(∇²Ibb)ᴰ*V*a
            bp.hessian[end,end] =bp.hessian[end,end] - dot(bp.r.*bp.bb.∇²I,bp.ϵ) # dot product
            copyto!(bp.x_hess_vec  , x)
        end
    end

    function evaluate_box_constraints(::BBPoint{N, Nx3, T} , 
                            emissivity_range::B=nothing, 
                            temperature_constraint::C = nothing) where {N, Nx3, T, 
                                                                        B <:Union{NTuple{2 , T} , Nothing} ,
                                                                        C <: Union{NTuple{2 , T} , Nothing} }
        
        return if isnothing(temperature_constraint) 
            DEFAULT_TEMPERATURE_RANGE[] 
        else 
            (temperature_constraint[1], temperature_constraint[2]) # limits on the BB temperature
        end
    end

    """
    fill_emissivity_box_constraint!(lb,ub,::VanderMatrix{N, CN, T, NxCN, CNxCN, P},
                val_bounds::NTuple{2,T}) where {N, CN, T, NxCN, CNxCN, P<:BernsteinSymPoly}

Evaluates box-boundaries for polynomial coefficients for `BernsteinSymPoly` 
polynomial basis
"""
    function fill_emissivity_box_constraint!(lb , ub , ::VanderMatrix{N, CN, T},
                    val_bounds::NTuple{2,T}) where {N, CN, T}

        fill!(lb , first(val_bounds))
        fill!(ub , last(val_bounds))

    end
    fill_emissivity_box_constraint!(lb , ub , V::VanderMatrix{N,CN,T} , ::Nothing) where {N , CN , T} = fill_emissivity_box_constraint!(lb , ub , V, T.(DEFAULT_EMISSIVITY_RANGE[]))

    function fill_emissivity_box_constraint!(lb , ub , ::VanderMatrix{N , CN , T},
                    val_bounds::NTuple{2 , <:Union{NTuple{CN , T}, StaticVector{CN , T}}}) where {N, CN, T}
                    
        copyto!(lb , first(val_bounds) ) 
        copyto!(ub , last(val_bounds) )
    end

    feval!(e::BBPoint , T::AbstractArray) = feval!(e,T[end])
    """
    feval!(e::BBPoint,t::Float64)

Evaluates bb intensity for temperature t
"""
function feval!(e::BBPoint , t::Number) # fills planck spectrum
        if t!=e.Tib[] # if current temperature is the same as the last recorded, 
            #a₁₂₃!(e_obj.amat,e_obj.λ,t) # filling amat
            Planck.a₁₂₃!(e.amat , e.λ , t) #fills amatrix
            Planck.ibb!(e.Ib, e.λ, e.amat) #fills BB spectrum
            e.Tib[] = t # save the current temperature
        end
        return e.Ib
    end
    residual!(e::BBPoint , T::AbstractArray) = residual!(e::BBPoint , T[end])
"""
    residual!(e::BBPoint,t::Float64)

Evaluates the residual vector between calculated and measured bb thermal 
emission intensity spectrum
Evaluates the residual vector between calculated and measured bb thermal 
emission intensity spectrum
"""
function residual!(e::BBPoint{N , N3 , T} , t::D) where {N , N3 , T <: Number , D <: Number}
        feval!(e , t)
        if t != e.Tri[] # if current temperature is the same as the last recorded, 
            e.ri .= e.I_measured .- e.Ib# calculating discrepancy
            e.r[] =0.5* norm(e.ri)^2 # discrepancy value
            e.Tri[]=T(t)# filling temperature of residual
        end
        return e.ri # returns residual vector
    end
    
    """
    disc(T , e::BBPoint)

Evaluates the least-square discrepancy between measured and calculates spectra
Input:
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature 
    e - (modified) current bb thermal emission point 
Evaluates the least-square discrepancy between measured and calculates spectra
Input:
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature 
    e - (modified) current bb thermal emission point 
"""
function  disc(T , e::BBPoint)
        residual!(e , T)# fills residuals
        return e.r[] # returns current value of discrepancy
    end

    ∇!(e::BBPoint , T::AbstractVector) = ∇!(e , T[end])

    """
    ∇!(e::BBPoint , t::Number)

Fills the first derivative of Planck function 
Input:
    T  - temperature 
    e - (modified) current bb thermal emission point   
Fills the first derivative of Planck function 
Input:
    T  - temperature 
    e - (modified) current bb thermal emission point   
"""
function ∇!(e::BBPoint , t::Number) # evaluates Planck function first derivative
        feval!(e,t)# refreshes amat and Ib
        if t != e.T∇ib[] # current temperature is not equal to the temperature of gradient calculation
            Planck.∇ₜibb!(e.∇I,t, e.amat,e.Ib)# fills Planck first derivative
            e.T∇ib[] = t # refresh gradient calculation temperature
        end
        return e.∇I
    end
    """
    grad!(g::AbstractVector,t::Float64 ,e::BBPoint)

In-place filling of the gradient of least-square problem of bb thermal emission spectrum fitting 
Input:
    g - gradient vector to be filled
    t  - temperature 
    e - (modified) current bb thermal emission ppoint 
"""
function grad!(g::AbstractVector, t::Number , e::BBPoint)
        g[] =  grad!(e , t)
        return g
    end
grad!(e::BBPoint , t::AbstractVector) = grad!(e , last(t))
grad!(g::AbstractVector, t::AbstractVector , e::BBPoint) = grad!(g, last(t) , e)
function grad!(e::BBPoint , t::Number)
        ∇!(e , t)
        residual!(e , t)
        if t != e.Tgrad[]
            ∇!(e , t)
            residual!(e , t)
        end
        return - dot(e.ri , e.∇I)
    end
    ∇²!(e::BBPoint , T::AbstractVector)=∇²!(e , T[end])

    """
    ∇²!(e::BBPoint , t::Number)

Fills the second derivative of Planck function 

    Input:
        T  - temperature 
        e - (modified) current bb thermal emission point  
"""
function ∇²!(e::BBPoint , t::Number)
        ∇!(e , t)# refreshes amat and Planck gradient
        if t != e.T∇²ib[]
           Planck.∇²ₜibb!(e.∇²I, t, e.amat, e.∇I) 
           e.T∇²ib[] = t # ref value
        end
        return e.∇²I
    end
    hess!(h , T::AbstractVector , e::BBPoint) = hess!(h , T[end],e)
    hess!(e::BBPoint , T::AbstractVector ) = hess!(e , T[end])
    """
    hess!(h,t::Float64,e::BBPoint)

In-place filling of least-square problem hessian matrix 
    Input:
        h - hessian 
        T - temperature 
        e - (modified) current bb thermal emission point  
        
In-place filling of least-square problem hessian matrix 
    Input:
        h - hessian 
        T - temperature 
        e - (modified) current bb thermal emission point  
"""
function hess!(h , t::Number , e::BBPoint{M , N , T}) where { M , N ,T <: Number} # calculates hessian of a simple Planck function fitting
        h[] = hess!(e , t)
        return h
    end
function hess!(e::BBPoint{M , N , T} , t::Number ) where { M , N ,T <: Number} # calculates hessian of a simple Planck function fitting
        if t != e.Thess[]
            e.Thess[] = T(t)
            ∇²!(e , t)
        end
        return dot(e.∇I , e.∇I) - dot(e.ri , e.∇²I)
    end



    function fitting_result(point::MWPPoint , results, probl , optimizer , ::Val{:full}) 
        return (T=temperature(point) , a=results.u[1:end-1],
                                    ϵ=point.vandermonde*results.u[1:end-1],
                                    res=results,
                                    problem = probl, 
                                    optimizer=optimizer)
    end
    fitting_result(point::Union{MWPPoint , BBPoint} , _, _ , _ , ::Val{:T}) = temperature(point)

    fitting_result(point::BBPoint , results , probl , optimizer ,::Val{:full}) = (T=temperature(point), 
                                                                    res=results, problem = probl ,  optimizer=optimizer)

    fitting_result(point::BBPoint , results , probl , optimizer::DefaultOptimizer ,::Val{:full}) = (T=temperature(point), 
                                                                    res=results, problem = probl ,  optimizer=optimizer)                                                                
    function trim_starting_vector_to_box!(v , lb , ub)
        for (i , (l , u)) in enumerate(zip(lb , ub))
             (l <= v[i]) && (v[i] <= u) ? continue :  v[i] = (l + u)/2
        end
    end

function fit_T!(p::Union{BBPoint , MWPPoint},
            o = DefaultOptimizer();
            emissivity_range::C=nothing, 
            temperature_range::B=nothing , 
            result_type::Val{D} = Val(:T)) where {B <: Union{AbstractVector , Nothing , NTuple{2}} , 
                                                  C <: Union{AbstractVector , Nothing , NTuple{2}} , 
                                                  D }  

        sv = get_default_starting_vector(p)
        return fit_T!(p , sv , o ;
                        emissivity_range = emissivity_range , 
                        temperature_range = temperature_range , 
                        result_type = result_type )
end
get_default_starting_vector(::BBPoint{M , N, T}) where {M , N, T} = MVector{1}(T(1000.0))
get_default_starting_vector(p::MWPPoint) = MVector(p.x)

function fit_T!(point::Union{BBPoint , MWPPoint}, 
                    starting_vector::AbstractVector,
                    optimizer = DefaultOptimizer();
                    emissivity_range::C=nothing, 
                    temperature_range::B=nothing , 
                    result_type::Val{D} = Val(:T)) where {B <: Union{AbstractVector , Nothing , NTuple{2}} , 
                                                  C <: Union{AbstractVector , Nothing , NTuple{2}} , 
                                                  D }
                                          
        (lb , ub) = evaluate_box_constraints(point, emissivity_range, temperature_range)
        trim_starting_vector_to_box!(starting_vector , lb , ub)
        (results , optimizer , problem) = _solve_problem(point , starting_vector , lb , ub , optimizer)
        return  fitting_result(point, results, optimizer , problem , result_type) 
                        
end

function (emp::Union{MWPPoint , BBPoint})(I::Union{AbstractVector , Number}; kwargs...) 
    set_measured!(emp , I)
    return emp(;  kwargs...)
end

function (emp::Union{BBPoint , MWPPoint})(; optimizer = DefaultOptimizer() , 
                                            starting_vector :: Union{Nothing , AbstractVector} = nothing,
                                            temperature_range = nothing , 
                                            emissivity_range = nothing , 
                                            result_type=Val(:T))

                    isnothing(starting_vector) && return fit_T!(emp , optimizer ; 
                                                            emissivity_range = emissivity_range , 
                                                            temperature_range = temperature_range , 
                                                            result_type = result_type)

                    return fit_T!(emp , starting_vector , optimizer ; 
                            emissivity_range = emissivity_range , 
                            temperature_range = temperature_range , 
                            result_type = result_type)
    end
"""
    try_emissivity(bb::BBPoint , ϵ::Number)

divides the measured intensity by the input emissivity 
"""
function try_emissivity(bbp::Union{BBPoint{N} , MWPPoint{N}}  , ϵ::ET  ) where ET <: Union{Number , SVector{N}} where N
    if ET <: Number
        (abs(ϵ) < 1e-8) && return 1e8
    else
        for e in ϵ
            (abs(e) < 1e-8) && return 1e8
        end
    end
    I_test = SVector(measured(bbp))
    # dividing by emissivity 
    T =  bbp( I_test./ϵ)
    set_measured!(bbp , I_test)
    return T
end
function _solve_problem(point , starting_vector , lb , ub , optimizer) error("To use multiwavelength pyrometry one must add Optimization package to the working env") end



function _solve_problem(point::MWPPoint{N, Nx3, P} , starting_vector , lb , ub , ::DefaultOptimizer) where {N, Nx3, P}
    
    return robust_lm_search!(point, 
                    MVector{P}(starting_vector) ,  
                    SVector(lb),  
                    SVector(ub))  

end

struct MultiWavelengthPyrometer{N , T , MWP} <: AbstractPyrometer{N,T}
    mwp::MWP
    function MultiWavelengthPyrometer(mwp::MWP) where MWP <: MWPPoint{N, Nx3, P, NxP, PxP, Pm1, NxPm1, Pm1xPm1, T} where {N, Nx3, P, NxP, PxP, Pm1, NxPm1, Pm1xPm1, T}
        return new{N , T , MWP}(mwp)
    end
end
MultiWavelengthPyrometer(λ::StaticVector{N , T}  ;  
                            xₒ::NTuple{P} = (0.5 , 0.5 , 0.5 , 1000.0) , 
                            i_measured::Union{Nothing , AbstractVector{T} , AbstractSpectralQuantity} = nothing,
                            polynomial_type = :bernsteinsym , 
                            I_sur::Union{StaticArray{Tuple{N}, T, 1}, Nothing} = nothing) where {N , P , T <: Number} = begin
    
    _i = if isnothing(i_measured) 
        MVector{N , T}(undef) 
    elseif isa(i_measured , AbstractVector) 
        MVector{N , T}(i_measured)
    elseif isa(i_measured , AbstractSpectralQuantity)
        MVector{N , T}(i_measured.(λ))
    end  
    
    mwp = MWPPoint( _i, 
                   MVector{N , T}(λ) , 
                   MVector{P , T}(xₒ) ; 
                   polynomial_type = polynomial_type , 
                   I_sur = I_sur  
                   )
    return MultiWavelengthPyrometer(mwp)                            
end
MultiWavelengthPyrometer{N}(l::AbstractVector;kwargs...) where N = MultiWavelengthPyrometer(MVector{N}(l); kwargs...) 

(p::MultiWavelengthPyrometer)(;kwargs...) = p.mwp(;kwargs...)
wavelengths(p::MultiWavelengthPyrometer) = p.mwp.bb.λ
measured(p::MultiWavelengthPyrometer) = measured(p.mwp)

calculated(p::MultiWavelengthPyrometer) = p.mwp.Ic
emissivity(p::MultiWavelengthPyrometer) = p.mwp.ϵ


measured(p::MWPPoint) = measured(p.bb)
measured(p::BBPoint) = p.I_measured

set_measured!(p::BBPoint , i::AbstractVector) = copyto!(p.I_measured , i)
set_measured!(p::MWPPoint , i::AbstractVector) = set_measured!(p.bb , i)
set_measured!(p::MultiWavelengthPyrometer , i::AbstractVector) = set_measured!(p.mwp , i)
"""
    lm_step(x, p::MWPPoint{N, Nx3, P}, λ::Real) where {N, Nx3, P}

Function to evaluate the Levenberg-Marquardt step 
"""
function lm_step(x, p::MWPPoint{N, Nx3, P}, λ::Real) where {N, Nx3, P}
    hess!(p, x) 
    H = SMatrix(p.hessian)
    J = SMatrix(p.jacobian)
    r = SVector(p.r)
    Jt = transpose(J)
    Δx = (H + λ * I) \ (- Jt* r)
    return Δx
end

function robust_lm_search!(mwp::MWPPoint{N, Nx3, P, NxP, PxP, Pm1, NxPm1, Pm1xPm1, T}, 
                x::MVector{P , T} ,  
                lower_bounds::SVector{P, T},  
                upper_bounds::SVector{P, T}; 
                atol::Real = 1e-8 , max_num_steps::Int = 500) where {N, Nx3, P, NxP, PxP, Pm1, NxPm1, Pm1xPm1, T}

    n = 0
    d_current = disc(x, mwp) 
    λ = 1e-2 
    x_trial = copy(x)
    
    while (d_current >= atol) && (n <= max_num_steps) 
        
        p_step = lm_step(x, mwp, λ) #lm step calculation 
        
        @. x_trial = x - p_step
    
        @. x_trial = clamp(x_trial, lower_bounds, upper_bounds)

        if any(isnan, x_trial)
            λ *= 10.0
            n += 1
            continue
        end 

        d_trial = disc(x_trial, mwp)

       if d_trial < d_current
            @. x = x_trial      # success
            d_current = d_trial 
            λ /= 5.0            # reduce demping factor 
        else
            λ *= 7.0            # unsuccesfull step increase demping
        end
        #@show λ
        n += 1
    end
    copyto!(mwp.x , x_trial)
    return (mwp , (λ = λ , n = n , d = d_current) , DefaultOptimizer())
end




"""
    fit_blackbody_safeguarded!(
    bb::BBPoint{N, Nx3, T_type}, 
    T_start::Real, 
    T_min::Real, 
    T_max::Real
            ) where {N, Nx3, T_type}


 Default solver for BBPoint with bc using Halley method
"""
function fit_blackbody_safeguarded!(
    bb::BBPoint{N, Nx3, T_type}, 
    T_start::Number, 
    T_min::Number, 
    T_max::Number ) where {N, Nx3, T_type}

    (T_min > T_max) && ((T_min , T_max) = (T_max , T_min))
    a = T_type(T_min)
    b = T_type(T_max)

    T_curr = clamp(T_type(T_start), a, b)
    
   # x_vec = MVector{1, T_type}(T_curr)
    max_iter = 30
    tol = T_type(1e-6)
    
    for iter in 1:max_iter

        f_val = grad!(bb , T_curr)  
        f_prime = hess!(bb, T_curr)                
        if abs(f_val) < tol
            break
        end
        
        if f_val > 0
            b = min(b, T_curr)
        else
            a = max(a, T_curr)
        end
        
        if (b - a) < tol
            T_curr = 0.5 * (a + b)
            break
        end
        
        f_prime_prime = 3.0 * dot(bb.∇I, bb.∇²I)
        
        # Halley: ΔT = (2 * f * f') / (2 * f'^2 - f * f'')
        denominator = 2.0 * f_prime^2 - f_val * f_prime_prime
        
        step_computed = false
        T_next = T_curr
        
        if abs(denominator) > 1e-12
            ΔT = (2.0 * f_val * f_prime) / denominator
            T_next = T_curr - ΔT
            
            if a + tol < T_next < b - tol
                step_computed = true
            end
        end
        
       # bisection if halleys doesnt work 
        if !step_computed
            T_next = 0.5 * (a + b)
        end
        
        T_curr = T_next
    end
    return (bb , nothing , DefaultOptimizer())
end


function _solve_problem(point::BBPoint , starting_vector , lb , ub , ::DefaultOptimizer)
    
    return fit_blackbody_safeguarded!(point, 
                    starting_vector[] ,  
                    lb[],  
                    ub[])  

end