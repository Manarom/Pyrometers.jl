"""
    stray_radiation_corrected_temperature(p::Pyrometer, Tmeasured::Number, source_intensity::Number)

Correct the measured temperature for a single-channel brightness pyrometer when the external 
background radiation (`source_intensity`) is provided directly as an integrated signal.

The baseline surface emissivity is automatically retrieved from the pyrometer instance. 
The system geometry is assumed to be an enclosure (`EnclosureGeometry`), meaning multiple mutual 
reflections are neglected, yielding a simple effective reflectivity of `r_eff = 1 - ϵ_surf`.
"""
function stray_radiation_corrected_temperature(p::Pyrometer, Tmeasured::Number, 
                                                    ϵ_surf::Number,
                                                    source_intensity::Number) 

        measured_signal = signal(p, Tmeasured , ϵ_surf)         
        r_eff   = one(ϵ_surf) - ϵ_surf
        reflected_signal = r_eff * source_intensity
        pure_signal = _extract_signals(p , measured_signal , reflected_signal)
        return p(pure_signal , ϵ_surf)
    end

"""
    stray_radiation_corrected_temperature(p::Pyrometer, Tmeasured::Number, 
                                                    ϵ_surf::AbstractSpectralQuantity,
                                                    source_intensity::Union{Number, IsothermalSpectralQuantity})

Surface emissivity can be temperature-dependent , external radiation must be provided as an isothermal quantity 
"""
function stray_radiation_corrected_temperature(p::Pyrometer, Tmeasured::Number, 
                                                    ϵ_surf::AbstractSpectralQuantity,
                                                    source_intensity::Union{Number, IsothermalSpectralQuantity}) 

        measured_signal = signal(p, Tmeasured , ϵ_surf)         
        r_eff   = SpectralReflectivity(ϵ_surf)
        reflected_signal = r_eff * source_intensity
        pure_signal = fix_temperature(_extract_signals(p , measured_signal , reflected_signal) , Tmeasured)
        return p(pure_signal , ϵ_surf) #
    end


"""
    stray_radiation_corrected_temperature(p::RatioPyrometer, Tmeasured::Number, 
                                                    source_intensity::NTuple{2 , T}) where T <: Number

Correct the measured temperature for a dual-channel ratio pyrometer when the external 
background radiation signals for both channels are directly provided as a tuple (`source_intensity`).

The baseline surface emissivities for both channels are automatically retrieved from the pyrometer instance. 
The system geometry is assumed to be an enclosure (`EnclosureGeometry`), meaning multiple mutual 
reflections are neglected, yielding independent channel reflectivities `r_eff = 1 - e_surf`. 
Returns the true temperature resolved from the corrected signal ratio.
"""
function stray_radiation_corrected_temperature(p::RatioPyrometer, Tmeasured::Number, e_surf::NTuple{2},
                                                    source_intensity::NTuple{2 , T}) where T <: Number

        e_surf1 , e_surf2 = e_surf
        si1 , si2 = source_intensity[1] , source_intensity[2]
        m_s1 , m_s2 = signal(p, Tmeasured , e_surf)         
        r_eff1 , r_eff2   = ( one(e_surf1) - e_surf1 , one(e_surf2) - e_surf2)
        r_s1 , r_s2 = r_eff1 * si1 , r_eff2 * si2 
        p_s1 , p_s2 = m_s1 - r_s1 , m_s2 - r_s2
        return p(p_s1/p_s2)
    end
"""
    stray_radiation_corrected_temperature(p::AbstractPyrometer, 
                                        Tmeasured::Number, 
                                        source_intensity
                                        ) where T<:Number

Takes emissivity from the pyrometer 
"""
stray_radiation_corrected_temperature(p::AbstractPyrometer, 
                                        Tmeasured::Number, 
                                        source_intensity
                                        )  = stray_radiation_corrected_temperature(p, Tmeasured,  
                                                get_emissivity(p) ,
                                                source_intensity
                                            )

"""
    stray_radiation_corrected_temperature(p::AbstractPyrometer, 
                                        Tmeasured::Number, 
                                        ϵ_surf::Union{T, NTuple{2, T}} ,
                                        source_intensity::AbstractContinuousOrDiscreteQuantity 
                                        ) where T<:Number

when the surface emissivity is a `Number` external radiation can be integrated 
"""
function stray_radiation_corrected_temperature(p::AbstractPyrometer, 
                                        Tmeasured::Number, 
                                        ϵ_surf::Union{T, NTuple{2, T}} ,
                                        source_intensity::Union{AbstractDiscreteQuantity , IsothermalSpectralQuantity};
                                        segbuf=nothing, rtol=sqrt(eps(Float64))
                                        ) where T<:Number 
                                        
    return stray_radiation_corrected_temperature(p, Tmeasured,  
                                                ϵ_surf ,
                                                integrate(p , Tmeasured , source_intensity , segbuf=segbuf, rtol=rtol)
                                        )
end
"""
    stray_radiation_corrected_temperature(p::AbstractPyrometer, 
        Tmeasured::Number, 
        ϵ_surf::Union{T , NTuple{2,T}},
        incident_radiation_function::AbstractSpectralQuantity ,  
        radiation_temperature::Number ; 
        segbuf=nothing, rtol=sqrt(eps(Float64)) ) where {T <: Number}

Correct the measured temperature when the external background source emits a known temperature-dependent 
intensity profile (which is not required to be an ideal blackbody spectrum).

"""
function stray_radiation_corrected_temperature(p::AbstractPyrometer, 
                                                Tmeasured::Number, 
                                                ϵ_surf::Union{T , NTuple{2,T}},
                                                incident_radiation_function::AbstractSpectralQuantity ,  
                                                radiation_temperature::Number ; 
                                                segbuf=nothing, rtol=sqrt(eps(Float64)) ) where {T <: Number}
        
        return stray_radiation_corrected_temperature(p , Tmeasured ,
                    ϵ_surf, 
                    integrate(p  , radiation_temperature ,  
                                incident_radiation_function , 
                                segbuf=segbuf, rtol=rtol)
        )
    end  




@inline function _extract_signals(::Pyrometer, measured, reflected)
    return measured - reflected
end

@inline function _extract_signals(::RatioPyrometer, measured::NTuple{2}, reflected::NTuple{2})
    m1, m2 = measured
    r1, r2 = reflected
    return (m1 - r1) / (m2 - r2)
end
    