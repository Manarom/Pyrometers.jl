   
#simple functions for intrepolation and numerical integration 

_local_interpolate(l::NTuple{N} ,  λ_nodes::AbstractVector, ϵ_nodes::AbstractVector) where N = ntuple(N) do i
        _local_interpolate(l[i] , λ_nodes , ϵ_nodes)
end
     """
    _local_interpolate(λ_target::Number, λ_nodes::AbstractVector, ϵ_nodes::AbstractVector)

Single point interpolation  of discrete data 
"""
@inline function _local_interpolate(λ_target::Number, λ_nodes::AbstractVector, ϵ_nodes::AbstractVector)
       
        idx = searchsortedfirst(λ_nodes, λ_target)
        
        if idx == 1
            return ϵ_nodes[begin]
        elseif idx > length(λ_nodes)
            return ϵ_nodes[end]
        end
        
        λ_start, λ_end = λ_nodes[idx-1], λ_nodes[idx]
        ϵ_start, ϵ_end = ϵ_nodes[idx-1], ϵ_nodes[idx]
        
        t = (λ_target - λ_start) / (λ_end - λ_start)
        return ϵ_start + t * (ϵ_end - ϵ_start)
    end
# internal function for trapz integration 
    """
    _trapz(x::AbstractVector, y::AbstractVector)

internal function for trapezoidal integration 
"""
function _trapz(x::AbstractVector, y::AbstractVector)
        N = length(x)
        @assert length(y) == N "Vectors must have the same length"
        N < 2 && return zero(eltype(y))

        s = zero(promote_type(eltype(x), eltype(y), Float64))
        
        @inbounds @simd  for i in 1:(N - 1)
            Δx = x[i+1] - x[i]
            ∑y = y[i+1] + y[i]
            s += Δx * ∑y
        end

        return s * 0.5
    end
    _simpson(x::AbstractRange , y::AbstractVector) =_simpson_even(x , y) 
    _simpson(x::AbstractVector , y::AbstractVector; is_evenly_spaced::Bool=false) = is_evenly_spaced ? _simpson_even(x , y) : _simpson_uneven(x , y)
    """
    _simpson_uneven(x::AbstractVector{TX}, y::AbstractVector{TY}) where {TX, TY}

internal function for Simpson's integration (uneven grids)
"""
function _simpson_uneven(x::AbstractVector{TX}, y::AbstractVector{TY}) where {TX, TY}
            N = length(x)
            @assert length(y) == N "Vectors must have the same length"
            N < 3 && return _tranpz(x , y)
            S = promote_type(TX, TY, Float64)
            s = zero(S)
            @inbounds @fastmath @simd for i in 1:2:(N - 2)
                h0 = x[i+1] - x[i]
                h1 = x[i+2] - x[i+1]
                h_sum = h0 + h1
                
                c0 = (2 * h0^2 + h0 * h1 - h1^2) / (6 * h0)
                c1 = (h_sum^3) / (6 * h0 * h1)
                c2 = (2 * h1^2 + h0 * h1 - h0^2) / (6 * h1)
                
                s += c0 * y[i] + c1 * y[i+1] + c2 * y[i+2]
            end
            if iseven(N)
                @inbounds @fastmath s += 0.5 * (x[N] - x[N-1]) * (y[N] + y[N-1])
            end

            return s
        end
        """
    _simpson_even(x::AbstractVector{TX}, y::AbstractVector{TY}) where {TX, TY}

Simpsons methods for evenly-spaced grids 
"""
    function _simpson_even(x::AbstractVector{TX}, y::AbstractVector{TY}) where {TX, TY}
        N = length(x)
        @assert length(y) == N "Vectors must have the same length"
        N < 3 && return _trapz(x, y) 

        S = promote_type(TX, TY, Float64)
        h = _get_step(x)

        last_simpson_idx = isodd(N) ? (N - 1) : (N - 2)

        s = S(y[begin]) + S(y[last_simpson_idx + 1])

        @inbounds @fastmath @simd for i in 2:last_simpson_idx
            coef = iseven(i) ? S(4) : S(2)
            s += coef * y[i]
        end

        integral_val = s * h / 3

        if iseven(N)
            @inbounds @fastmath integral_val += 0.5 * h * (y[N-1] + y[N])
        end
        return integral_val
    end
    _get_step(x::AbstractRange) = step(x)
    _get_step(x::AbstractVector)  =  x[2] - x[1]