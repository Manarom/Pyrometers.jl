### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 5e712312-0fc7-4205-84cc-834d57b814a3
begin 
	import Pkg 
	notebook_dir = @__DIR__()
	Pkg.activate(notebook_dir)
	Pkg.resolve()
	using Revise
	using Pyrometers  , Plots , PlutoUI , PrettyTables , DelimitedFiles , Interpolations 
	using QuadGK , Pyrometers.StaticArrays
	src_dir = joinpath(abspath(joinpath(notebook_dir,"..")),"src")
end;

# ╔═╡ 0bee71f4-5961-4c80-8592-2b6c0d1b58a8
	begin 
		using ForwardDiff
		using Optim , NLSolversBase
	end

# ╔═╡ 05c05c84-02d4-4b7f-83df-bd1fa3e4ee4d
using BenchmarkTools , Test

# ╔═╡ a8a27a41-95a2-41eb-a1cf-d3d51b2ec52e
using LinearAlgebra

# ╔═╡ 30743a02-c643-4bdc-837e-b97299f9520a
md"""
#  `Pyrometers.jl` package usage  

##### Package **`Pyrometers.jl`** contains methods for virtual radiation pyrometers of different types. 
____________________
### Installation

To run this notebook, you need:
1) Install `julia` language itself from its official [download page](https://julialang.org/downloads) 
2) Install [Pluto](https://plutojl.org/) notebook from `julia` REPL by entering the following commands line-by-line:
```julia
import Pkg
Pkg.add("Pluto")
using Pluto
Pluto.run()
```
The last line will launch the Pluto starting page in your default browser 
3) Copy the entire GitHub [repository](https://github.com/Manarom/Pyrometers.jl.git) to your local folder
4) Open this notebook file located at  `project_folder\notebooks\radiation_pyrometry.jl` in `Pluto` by providing the full path to the *"Open a notebook"* text field on `Pluto`'s starting page.

"""

# ╔═╡ f728a59d-c78c-45af-a9e5-656be490eb4f
Revise.retry()

# ╔═╡ abdc809b-b53c-4dff-ba6f-c636c73f3fca
const Planck = Pyrometers.Planck

# ╔═╡ bf833e74-f9e7-4b60-b6bc-2a6a58c5c901
const MY_GLOBAL_SEGBUF = QuadGK.alloc_segbuf(Float64, Float64, Float64, size=5000);

# ╔═╡ 171409eb-22b5-4bc5-a8e2-eac0932a24f3
PlutoUI.TableOfContents(indent=true, depth=4, aside=true)

# ╔═╡ 643d9ff3-3a09-46c9-9013-92d111ccb229
plot_common_args = (grid = true, gridlinewidth=3, gridstyle = :dot,minorgrid=true, box = :on, linewidth = 3);

# ╔═╡ d5ee3913-66be-47d7-a755-699ba64b4f98
md"""
### Introduction

This notebook demonstrates examples of using the **Pyrometers.jl** package. 

The purpose of this package is to simulate virtual pyrometers of **five different types**:

* `SpectralBandPyrometer`: A classical pyrometer operating within a fixed wide spectral band.
* `TwoBandsRatioPyrometer`: A spectral ratio pyrometer that accounts for the finite width of two channels.
* `SingleWavelengthPyrometer`: A standard pyrometer operating at a single, fixed wavelength.
* `TwoWavelengthRatioPyrometer`: A classical spectral ratio pyrometer that does not account for the finite width of its spectral channels.
* `MultiwavelengthPyrometer`: A pyrometer based on least-squares fitting of discrete spectral intensity data across multiple wavelengths (theoretically eliminating the need to know surface emissivity).

All of these pyrometers implement **five core functions**:

* `measure`: Converts input intensity into temperature.
* `signal`: Converts temperature into a signal (the inverse of `measure`).
* `convert_temperature`: Transforms one temperature to another (e.g., to account for true surface emissivity).
* `stray_radiation_corrected_temperature`: Excludes external radiation specified as an intensity.
* `external_source_corrected_temperature`: Excludes external radiation specified as a source with a given temperature and emissivity.

All these methods are designed to work with simple scalar values (numbers), but they can also account for wavelength- and temperature-dependent surface emissivity, as well as external radiation sources and intensities. 

To handle spectral properties, **Pyrometers.jl** provides a set of types that extend the `AbstractSpectralQuantity` types from the **PlanckFunction.jl** package. It also implements a simple "algebra-on-spectra" approach to perform basic mathematical operations on properties that depend on wavelength and (optionally) temperature.

"""

# ╔═╡ d442014a-20e6-4be4-ac7f-f13de329dec5
md"""
### I. `Blackbody` vs `real surface` thermal emission 
_______________________

All heated bodies emit thermal radiation. According to Planck's law, the spectrum of ideal emitter (called the *blackbody*) is governed solaly by its temperature. The blackbody spectral intensity can be calculated as follows \


``I_{blackbody}(\lambda , T) =  \frac{C_1}{\lambda ^5} \cdot \frac{1}{e^{\frac{C_2}{\lambda T } } - 1}``, \


where ``C_1`` = $(Planck.C₁), ``W \cdot μm/m² \cdot sr`` and ``C_2`` = $(Planck.C₂), ``μm \cdot K`` , ``\lambda`` - wavelength in ``\mu m``, ``T`` - temperature in Kelvins

A real surface thermal emission intensity is lower than the one of the blackbody. The fraction of blackbody thermal radiation intensity emitted by a real surface is characterized by directional spectral emissivity ``\epsilon (\lambda, T,\vec{\Omega})``:

``I_{real\ \ surface}(\lambda , T, \vec{\Omega} ) = \epsilon (\lambda, T,\vec{\Omega}) \cdot  I_{blackbody}(\lambda , T)``, \


here  ``\vec{\Omega}`` stays for direction. \


It is interesting that, unlike the blackbody, the real surface thermal emission (in general) depends  on the direction of radiation. Therefore, the most general characteristic for thermal radiation of a real surface is the `directional spectral emissivity`.

In [PlanckFunctions.jl](https://manarom.github.io/PlanckFunctions.jl) package there are several functions to calculate the blackbody thermal emission spectra (and various derivatives, integrals etc.).
The following figure show the impact of spectral emissivity on the real surface thermal emission intensity.
"""

# ╔═╡ 27b3c586-9eb0-4a51-b9ca-a9c0379fccdf
@bind  λ_BB PlutoUI.combine() do Child
	md"""
	Blackbody spectral range, ``\mu m`` : \
	``\lambda_{left}`` = $(
		Child(Slider(0.1:0.1:50,default=0.1,show_value = true))
	)   -- 
	 $(
		Child(Slider(0.1:0.1:50,default=15.0,show_value = true))
	)  ``\lambda_{right}`` 
	"""
end

# ╔═╡ f22d22b6-5d98-4cc4-998f-a53e92809618
@bind  T_BB PlutoUI.combine() do Child
	md"""
	Blackbody temperatures: \
	T₁ = $(
		Child(Slider(-273.0:1:2500,default=800,show_value = true))
	) ``^o C`` \
	T₂ = $(
		Child(Slider(-273.0:1:2500,default=900,show_value = true))
	)  ``^o C`` 
	"""
end

# ╔═╡ 7cc110e9-7655-4dfc-b1e0-ab3905866425
@bind  scales_BB PlutoUI.combine() do Child
	md"""
	Scales : \
	xscale = $(
		Child(Select([:identity,:ln,:log10]))
	)   yscale  
	 $(
		Child(Select([:identity,:ln,:log10]))
	)   
	"""
end

# ╔═╡ 59b01938-9235-4d39-be4e-e4fdf9f94418
md""" 
### II. *Algebra* of spectra
_______________________

In real life the emissivity is both wavelength and temperature dependent (in general it also depends on direction).  **Pyrometers.jl** provides special types for spectral properties (they all are subtypes of `AbstractSpectralQuanity`): 

"""

# ╔═╡ 6e3f01fb-a0bb-438b-a95a-2c503cc4961e
Markdown.parse(join(["$(i). $(last(split(string(t), ".")))" for (i, t) in enumerate(subtypes(Pyrometers.AbstractSpectralQuantity))], "\n"))

# ╔═╡ b7117f7b-ad22-4b09-8aa3-3f6a3878c149
md""" 

**AnalyticalSpectralQuantity** is a wrapper around spectral quantity which has analytical expressions for the quantity itself and its first and second derivatives with respect to temperature:

```julia 
e_a = AnalyticalSpectralQuantity(f , df , d2f)
# each of f , df , d2f is a function of two arguments l , t, wavelength and temperature respectively
e_a(l , t) # returns function value
eval_Dₜ(e_a , l , t) # returns the tuple of (value , df , d2f)
```

**IsothermalSpectralQuantity** wrapper for temperature independent quantities:

```julia 
e_i = IsothermalSpectralQuantity(f ) # f - function of wavelength in μm
e_i(l)
e_i(l , t) # returns function value (ignores the second argument)
eval_Dₜ(e_i , l , t) # returns the tuple of (value , 0.0 , 0.0)
```
**GenericDifferentiableSpectralQuantity** is a wrapper around function which can be differentiated by second argument using `ForwadDiff`

"""

# ╔═╡ 5efd8d48-7fff-48fd-b361-2f06cd8bad53
begin  
	e_analyt = Pyrometers.AnalyticalSpectralQuantity((l,t)-> 2.0e-4*t + 5e-3*l , 
													(l,t)-> 2.0e-4 , 
													(l,t)-> 0.0)
	e_iso = Pyrometers.IsothermalSpectralQuantity(l-> 2.0e-1 + 1e-3*l)
	e_dif = Pyrometers.GenericDifferentiableSpectralQuantity((l,t)-> 2.0e-4*t + 5e-3*l)
end;

# ╔═╡ af619a4d-b700-4996-960c-c9d0e75eac6a
md"""
Select quantity to show  $(@bind quant_to_show Select([1=>"f" ,  2=>"df" , 3=> "d2f"]))
"""

# ╔═╡ 0aad4eb6-e5dc-40e3-82b6-4d64f5b1af4f
md""" The following figure shows spectra for selected quantities"""

# ╔═╡ e9e16b2b-88d7-40b6-997f-3a40e673faa8
begin 
	p_quant = Plots.plot(;plot_common_args...)
	
	for t = (1000.0 , 1200.0)
		l = range(1.0 , 2.0 , 100)
		for (q , n) in zip((e_analyt , e_iso , e_dif) , ("Analytic" , "Isothermal" , "Differenciable"))
			f = Base.Fix1(Pyrometers.Planck.eval_Dₜ , q)
			o = [ f(l_i , t)[quant_to_show] for l_i in l]   
			plot!(p_quant , l ,o;plot_common_args... , label = n*": T = $t")
		end	
	end
	p_quant
end

# ╔═╡ 63d5b5bb-5e7d-4e83-83ae-3972f2465e28
md"""
**PlanckEmitter** is a Blackbody source wrapped into `AbstractSpectralQuanity` subtype

```julia 
bb_source = PlanckEmitter()
bb_source(l , T) # returns spectral intensity for bb source with temperature t at wavelength l
```
"""

# ╔═╡ 95baff6c-2af3-401d-9db3-d6e6e31ab70a
md""" Other quantities are used to implement simple algebra, by wrapping two `AbstractSpectralQuantities` into a single object ,  mathematical operators `*` , `/` , `+` when applied on two `AbstractSpectralQuantity` objects or any `Number` and `AbstractSpectralQuantity` combination return the object of **SpectralQuantitiesProduct** , **SpectralQuantitiesRatio** , **SpectralQuantitiesSum** respectively 

```julia

bb = PlanckEmitter()
ϵ = GenericDifferentiableSpectralQuantity((l,t)-> 2.0e-4*t + 5e-3*l)

i = ϵ * bb # returns SpectralQuantitiesProduct
i(l , t) # returns ϵ(l , t) * bb(l , t)
eval_Dₜ(i, l , t) # returns ϵ(l , t) * bb(l , t) values and first and second derivatives 
```

`AbstractSpectralQuantity`'s can also be combined with numbers:

```julia

e_new = 1 - 2 * ϵ # returns callable spectral quantity

```

This algebraic constructions can be usefull together with PlanckEmitter by imitating real surfaces emission , e.g.:

```julia
bb = PlanckEmitter()
e1 = GenericDifferentiableSpectralQuantity(f) # here f is any autodifferenciable function 
e2 = GenericDifferentiableSpectralQuantity(f2) # here f2 is any autodifferenciable function 
i = bb * (e1*e2)/(e1 + e2) # returns an object of `AbstractSpectralQuantity` subtype 
```

"""

# ╔═╡ 27203e6b-9283-4da3-beb2-9f671100860f
md" Adjust temperature $(@bind _t Slider(200:1:1700 , show_value = true , default = 1000))"

# ╔═╡ 5010ec7e-b83c-47c6-8c84-17bbcb8e40a1
begin 
	bb_source = PlanckEmitter()
	e_test  = GenericDifferentiableSpectralQuantity((l,t)-> 1.0e-4*t + 5e-3*l)
	i_combined =  2*(1 - e_test)*bb_source
	p_combined = Plots.plot(;plot_common_args...)
	
	for t = (_t  , )
		l = range(0.5 , 5.0 , 100)
		f = Base.Fix1(Pyrometers.Planck.eval_Dₜ , i_combined)
		v = [ f(l_i , t)[1] for l_i in l]   
		dv = [ f(l_i , t)[2] for l_i in l]  
		d2v = [ f(l_i , t)[3] for l_i in l]  
		plot!(p_combined , l , v;plot_common_args...  , label = "i")
		plot!(p_combined , l , dv;plot_common_args...  , label = "didt")
		plot!(p_combined , l , d2v;plot_common_args...  , label = "d2idt2", yscale = :log10)
	end
	p_combined 
end

# ╔═╡ 8a066ee5-80e9-462f-9a61-15851468aa63
md"""
### III. Partial radiation pyrometry
_______________________

As far as the `blackbody` thermal radiation energy strongly depends on temperature, this quantity can be used to measure the temperature of a real surface. This is the general idea of partial radiation pyrometry: **measure intensity to get the temperature**. As far as the intensity is a directional quantity, a pyrometer needs collimating optics (a telescope).The real surfaces emissivity often varies sufficiently with the wavelength, at the same time, partial radiation pyrometers assume constant emissivity (so-called `grey`-band approximation). Thus, for industrial purposes, it is useful to have several pyrometers, each working within a relatively narrow spectral band. In the spectral range of a partial radiation pyrometer, emissivity should not vary significantly to make the assumption of constant emissivity relevant.

"""

# ╔═╡ b7fac177-c211-4635-992f-e6473be7bdae
md"""
Dictionary **Pyrometers.DefaultPyrometersTypes** contains default pyrometers type names together with spectral range. Custom pyrometer can be created by providing its type (name), wavelength or wavelentgh range, working emissivity: 
``` julia
# creating custom pyrometer objects
p_custom = RadiationPyrometers.Pyrometer(type = "Custom",λ = [2.5, 3.7],ϵ=0.65)
```
"""

# ╔═╡ e2a9aa39-2490-4681-89d3-a01f058f6feb
pretty_table(HTML,Pyrometers.DefaultPyrometersTypes,top_left_string ="Table of default pyrometers types provied by `RadiationPyrometers.jl` package and corresponding wavelength regions",wrap_table_in_div=true)


# ╔═╡ 02eee968-ff43-4b82-8d68-efede1a220dd
md"""
After creating the **Pyrometer** object, it can be used to "measure" the temperature from signal, according to the pyrometer's spectral range using **Pyrometers.measure(pyr,measured_intensity)** function (the intensity should be provided in correct units [W/m²⋅sr⋅μm]) or by calling pyrometer directly (all pyrometers are functors). 

```julia
	p = SpectralBandPyrometer(2.0 , 3.0)
	i = fix_temperature(PlanckEmitter() , 1234.6)
	p(i) #returns  1234.6
```

"""

# ╔═╡ e9d0e216-07b5-4e41-925a-0c530e6ddf7b
md"""
 In practice it is quite common situation when one needs to find the emissivity in pyrometer range from real temperature of the surface at some point and temperature measured by the pyrometer, in this case the emissivity of virtual pyrometer can be adjusted to make the "measured" temperature be equal to the real one. This can be done using.

```julia
	Pyrometers.fit_ϵ!(p::Pyrometer,Tmeasured::Float64,Treal::Float64) 
```

This function adjusts the emissivity of pyrometer object.

"""

# ╔═╡ 9b08b767-7e8f-4483-9f2f-226022ce10e4
md"""
	If emissivity of partial radiation pyrometer is incorrect,  the measured temperature is also inaccurate. It is interesting to look how various pyrometers (with their emissivity set to one) `measure` the temperature of a real surface. The following figure shows the real surface thermal emission intensity and several common pyrometers types working regions. In the legend their `measured` temperature is shown. The `mesured` temperature for each pyrometer type is obtained by fitting the blackbody power to the real surface power both integrated over pyrometer's working spectral range.  
	"""

# ╔═╡ 144b40ea-71c7-421f-8117-eab267ea5daf
Pyrometers.produce_pyrometers()

# ╔═╡ 712828a7-fb54-42e6-95fc-233243190f59
md"Real surface temperature $(@bind T_pyr Slider(range(10,3000,1000),default=1500,show_value=true) ) "

# ╔═╡ f763d449-2a7a-4008-a183-823a774bc25e
#savefig(plot_pyrometers,joinpath(notebook_dir,"Pyrometers.png"));

# ╔═╡ 667f7c30-56e0-461f-b35b-c924007eb9f2
md"""
For this particular material the type -`F` pyrometer readings are closer to the real temperature, because of the real surface emissivity being closer to one for this pyrometer's spectral region. All results are summarized in the following table.
"""

# ╔═╡ d08ec8f2-7689-4043-9f28-da06ab0124b9
md"""
### IV. Blackbody reference source emissivity
_______________________

A real-life pyrometer does not directly measure radiation intensity due to the spectral dependence of its sensitivity and various other factors; therefore, it requires calibration. To calibrate the pyrometer, it is placed in front of a reference source that closely approximates blackbody thermal emission. Typically, this blackbody reference is a specially designed furnace with a heated cavity and a precise temperature controller.

During the calibration process, the pyrometer operator measures the temperature of the reference source. Ideally, the measured temperature should match the temperature set on the reference controller; however, in practice, these values rarely coincide exactly. This discrepancy arises from the non-ideal nature of the reference source. The spectral emissivity of a real reference is not exactly equal to one and varies with both wavelength and temperature. To account for this deviation, the reference source is supplied along with a calibration table, which looks like the following:

"""

# ╔═╡ 0c9fe7b1-374c-4fb8-9cfe-9337389713bf
md"""
The first row in the table represents the reference temperatures, while the subsequent columns show temperatures measured by different pyrometers. Each column is labeled according to the pyrometer type. It should be noted that the temperatures listed in the table vary both with the pyrometer type (across each row) and the temperature values themselves (down each column). This indicates that the spectral emissivity of the reference source depends on both the wavelength (corresponding to the pyrometer type) and the temperature. Since the various pyrometer types collectively cover a broad spectral range from 2 to 14 μm, the measured temperature data can be utilized to determine the spectral emissivity of the reference and its temperature dependence from the calibrations temperatures table provided above.

The following figure shows the spectral emissivity of the blackbody reference, calculated from the temperature calibration table shown above.
"""


# ╔═╡ 72947d97-0a97-4064-a2fe-08d19dec0f0e
md"""
 ### V Two temperature emission
"""

# ╔═╡ 9ce196c1-8915-46da-9aba-f13d7655959a
md"""
 	Now the emission is proportional to the spectral intensity of two mixed BB sources with different temperatures.

``I_{meas} = \epsilon\mathcal{b}(\lambda , T_1) + (1 - \epsilon)\mathcal{b}(\lambda , T_2)``

Let's check if it is possible to measure the surface temperature in a presence of external emission with a much higher temperature
"""

# ╔═╡ 3581aa29-714b-422a-8feb-d1a0c3ebeec7
begin 
	data_folder = joinpath(notebook_dir , "data")
	data_names = readdir(data_folder)
end;

# ╔═╡ b9bee300-59a4-4c7a-b525-439f5c62253e
md""" Select material: $(@bind material_selection Select(data_names , default = "zrb2.txt"))"""

# ╔═╡ 2ad3ec82-54a2-49ac-94ef-579f808dfb1a
begin 
	rt_emissivity_data = readdlm(joinpath(data_folder , material_selection)) # loading file, actually this file is two-column data with no headers, but this is ok
	rt_emissivity_interpolation = linear_interpolation(rt_emissivity_data[:,1],rt_emissivity_data[:,2],extrapolation_bc = Interpolations.Flat())
end;

# ╔═╡ 4d6337aa-cfc7-4154-a395-5aa53e23d01a
begin 
	T₁ = T_BB[1] + Planck.Tₖ # converting to Kelvins
	T₂ = T_BB[2] + Planck.Tₖ # converting to Kelvins
	λ_bb = collect(range(λ_BB...,length=500));
	ibb1 = Planck.ibb.(λ_bb,T₁ );ibb2 =  Planck.ibb.(λ_bb,T₂)
	
	p_bb = Plots.plot(λ_bb,ibb1,label="blackbody T=$(T₁),K", xscale=scales_BB[1],yscale =scales_BB[2],fillrange=0, fillalpha=0.3)
	
	plot!(λ_bb,ibb2,label="blackbody T=$(T₂),K", xscale=scales_BB[1],yscale =scales_BB[2],fillrange=0, fillalpha=0.3)
	
	Plots.plot!(λ_bb,ibb1.*rt_emissivity_interpolation(λ_bb),label="real surface T=$(T₁),K", xscale=scales_BB[1],yscale =scales_BB[2],fillrange=0, fillalpha=0.2)	
	
	Plots.plot!(λ_bb,ibb2.*rt_emissivity_interpolation(λ_bb),label="real surface T=$(T₂),K", xscale=scales_BB[1],yscale =scales_BB[2],fillrange=0, fillalpha=0.2)
	xlabel!("Wavelength, μm")
	ylabel!("Spectral intensity, W/m²⋅sr⋅μm")
end

# ╔═╡ 03d76e64-ebf4-432b-b9be-d4cb26275f55
begin 
	e_real_BB = rt_emissivity_interpolation(λ_bb) # this spectral emissivity was measured up to 18 μm, thus for higher wavelengths it uses flat extrapolation
	Plots.plot(λ_bb, e_real_BB, label=nothing, linewidth=3.0)
	title!("Real (experimental) surface spectral emissivity")
	xlabel!("Wavelength, μm");ylabel!("Spectral emissivity (ϵ)")
end

# ╔═╡ c69acbf6-94fb-4ac3-8d56-d1f9dda11440
begin
	λ_pyr = collect(range(0.1,18,1000))
	pyrometers_vector = sort(Pyrometers.produce_pyrometers())# returns a vector of all default pyrometers
	N = length(pyrometers_vector) + 1
	# calculationg the real surface thremal radiation spectrum
	
	real_i = Planck.ibb.(λ_pyr,T_pyr).*rt_emissivity_interpolation(λ_pyr)
	
	data_legend = Matrix{String}(undef,N,1)
	data_legend[1] = "T_real =$(round(T_pyr))"
	# poltting surface thremal emission spectrum
	
	plot_pyrometers = Plots.plot(λ_pyr,real_i,xscale=scales_BB[1],yscale =scales_BB[2],fillrange=0, fillalpha=0.3,dpi=600,label = data_legend[1],legend_background_color=:white,legend_foreground_color = :black,legend_position=:right)

	real_i_interp = linear_interpolation(λ_pyr,real_i)
	xlabel!("Wavelength , μm")
	ylabel!("Thermal radiation intensity")
	#ylabel!()
	max_val = maximum(real_i)
	t_em_unity = Vector{Float64}(undef,length(pyrometers_vector))
	t_em_acttual = Vector{Float64}(undef,length(pyrometers_vector))
	
	for j in 1:length(pyrometers_vector)# iterating over virtual pyrometers vector
		# the following function checks if current pyrometer is narrow band
		
		ppp = pyrometers_vector[j]
		is_two_wavelength_pyrometer = Pyrometers.is_spectral_band(ppp)

		l_cur = is_two_wavelength_pyrometer ? ppp.λ : [ppp.λ[]-0.2,ppp.λ[]+0.2 ]
		
		λ_pyr_interp = collect(range(l_cur...,length=30))
		# calculating the measured by the pyrometer value 
		measure_intensity =  is_two_wavelength_pyrometer ? Pyrometers._simpson(λ_pyr_interp,real_i_interp(λ_pyr_interp)) : real_i_interp(ppp.λ[1])
		
		# setting emissivity to one
		Pyrometers.set_emissivity!(ppp , 1.0)
	
		measured_temp =ppp(measure_intensity)
		# temperature measured by the current pyrometer 
		data_legend[j+1] = "$(ppp.type) : T="*string(round(measured_temp))
		# plotting current pyrometer spectral range
		region_flag = 
		plot!(l_cur,[max_val,max_val],fillrange=0, fillalpha=0.5,label=data_legend[j+1])
		# remember the value of temperature with unit emissivity
		t_em_unity[j] = measured_temp 
		# calculating the averaged gray-band emissivity
		Pyrometers.fit_ϵ!(ppp,measured_temp,T_pyr)
		# calcaulting the averaged emissivity wthin the pyrometers spectral band (or at fixed wavelength)
		 t_em_acttual[j] =  round(ppp(measure_intensity))
		 
	end
	plot!(twinx(),λ_pyr,rt_emissivity_interpolation(λ_pyr),linewidth=4,linecolor=:red,label=nothing,alpha=0.3,ylabel ="Real surface spectral emissivity" )
	
	
	# emissivities table 
	data = Matrix{Any}(undef,length(pyrometers_vector),5)
	e_grey = [p.ϵ[] for p in pyrometers_vector]
 	data[:,3:end] .= hcat(t_em_unity,e_grey, t_em_acttual)
	data[:,1] .= [p.type for p in pyrometers_vector]
	data[:,2] .= [string(p.λ) for p in pyrometers_vector]

end;

# ╔═╡ a861d56f-f6c9-4754-b7a9-ed63713f1f2f
	pyr_table = pretty_table(HTML,data, column_labels= ["type","λ region,μm","T₀ (ϵ=1),K","grey-ϵ", "T₁ (grey-ϵ),K"],top_left_string ="Table of temperatures `measured` by different pyrometers  with the spectral emissivity settled to one (T₀) and to the calculated grey-ϵ and temperature measured after setting gray band emissivity to the right value (T₁) the real temperature is Tᵣ=$(T_pyr)" )

# ╔═╡ 7071a6f4-e296-4e53-8e6f-24f1f038c1a5
plot_pyrometers

# ╔═╡ c5ac80ee-8143-4c28-bbff-2ac761c71fac
begin
	bb_calibration_table_data = readdlm("BBethalon")
	table_header =["Tref";]
	all_types = [getfield(p,:type) for p in pyrometers_vector]
	table_header = vcat(table_header,all_types)
	pr_tbl = pretty_table(HTML,bb_calibration_table_data,column_labels= table_header,top_left_string = "Example of table data for the blackbody reference source (all tempeatures are in Celsius)")
	bb_calibration_table_data .+= Planck.Tₖ # converting table data to 
	ref_T = @view bb_calibration_table_data[:,1] # reference source temperature
	pr_tbl
end

# ╔═╡ bc2d93ae-6c30-462c-96e0-30fdb84d7c63
md"""
Adjust blackbody reference temperature 

$(@bind T_ref1  Slider(ref_T,show_value=true,default = ref_T[1])) 

$(@bind T_ref2  Slider(ref_T,show_value=true,default = ref_T[1])) 

$(@bind T_ref3  Slider(ref_T,show_value=true,default = ref_T[end])) 

"""

# ╔═╡ d3199b6e-9779-4def-b701-fe85d1035045
begin 
	full_wavelengths_range = Pyrometers.full_wavelength_range(pyrometers_vector)
	jj = indexin(T_ref1,ref_T)[]
	foreach(pyrometers_vector) do p
		Pyrometers.set_emissivity!(p , 1.0)
	end

	
	ej = Pyrometers.fit_ϵ_wavelength!(pyrometers_vector,T_ref1,bb_calibration_table_data[jj,2:end])
	pppp = Plots.plot(full_wavelengths_range, ej, title="Spectral emissivity of the blackbody reference",label ="T = $(ref_T[jj])",grid=true)
	for T_reference in (T_ref2 , T_ref3)
		foreach(pyrometers_vector) do p
			Pyrometers.set_emissivity!(p , 1.0)
		end
		global jj = indexin(T_reference,ref_T)[]
		global ej = Pyrometers.fit_ϵ_wavelength!(pyrometers_vector , ref_T[jj] , bb_calibration_table_data[jj , 2:end])
		Plots.plot!(pppp , full_wavelengths_range, ej, label=" T = $(ref_T[jj])")
	end
	xlabel!(pppp,"Wavelength, μm")
	ylabel!(pppp,"Emissivity")
	pppp
end

# ╔═╡ 8a6fe87d-0f7c-4577-9ac2-ea1ecf71016b
pyrometers_vector

# ╔═╡ 86af6afc-b28a-4e84-952a-bd29710374f8
λ2 = collect(range(0.1,15.0,1000));

# ╔═╡ f181980f-bf72-4468-8daa-9461c6c901e0
md"""
Select the material (or fixed emissivity) : $(@bind  emissivity_type Select(vcat(["fixed"] , data_names), default = "fixed"))

"""

# ╔═╡ 36ba2396-bb5e-4d22-a58e-9ab27cd18b2d
md" Sample emissivity , ϵ = $(@bind ϵ_fixed  Slider(1e-4:1e-4:1.0 , show_value = true , default = 0.5))"

# ╔═╡ efc35420-d0e6-4795-94b6-d43289b4de44
begin 
	if emissivity_type == "fixed"
		eint = Returns(ϵ_fixed)
	else
		_data = readdlm(joinpath(data_folder , emissivity_type))
		#=if contains(emissivity_type , "quarz")
			@. _data[:,2] = 1.0 - _data[:,2]
		end=#
		eint = @views linear_interpolation(_data[:,1] , _data[:,2] , extrapolation_bc=Line())
	end
	e_int_iso = Pyrometers.IsothermalSpectralQuantity(eint)
end;

# ╔═╡ 91bbd553-4e4a-431d-9d54-b0f4882fd426
md"""
Sample temperature  ``T_1 `` = $(@bind T1 Slider(300.0:1e-2:3000 , show_value = true , default = 1000.0)), K

"""

# ╔═╡ 15f1519b-d924-4fa9-b212-eebba75c544a
md"""
External radiation temperature ``T_2`` = $(@bind T2 Slider(300.0:1e-2:4000 , show_value = true , default = 2800.0)), K
"""

# ╔═╡ 8cef05a1-2974-4c38-b73e-fa706f347fcc
md" Show wavelength range: $(@bind λ_show RangeSlider(range(extrema(λ2)... , 1000)))"

# ╔═╡ 1811e43c-f7db-47b1-9b83-bb38455d7db3
pyrometers_vector2 = deepcopy(pyrometers_vector);

# ╔═╡ 54339700-71fd-48bf-a2ef-0c3267b9d81b
md" **Recalculate T matrix $(@bind is_recalculate CheckBox(false))**"

# ╔═╡ 7f76bc22-77c3-4eb3-9d49-588e653df2e7
md"""

###### The following section calculates the `measured` temperature and relative error of temperature mesurement for selected or custom pyrometer type   

"""

# ╔═╡ a0197a9a-34bf-4a3e-af8a-c23ea777f482
md" Use custom pyrometer : $(@bind use_custom CheckBox(default = false))"

# ╔═╡ 7793976e-6714-4bc1-9d97-412dd2a67480
@bind custom_waves PlutoUI.combine() do Child

md"""
	
Custom pyrometer wavelength range:

λleft = $(Child("left", NumberField(0.1:1e-3:20 , default = 8.0)))

λright = $(Child("right", NumberField(0.1:1e-3:20 , default = 9.7)))

"""
end

# ╔═╡ e480137d-b6d9-4e18-92f0-640292bbb5f0
function two_planck(l , ϵ , T1 , T2)
	return ϵ * Planck.ibb(l , T1) + (1.0 - ϵ) * Planck.ibb(l , T2)
end

# ╔═╡ 6342e92b-4434-4e4b-aa2f-56405277caed
begin 
	ϵ  = e_int_iso.(λ2)
	I1 = @. ϵ * Planck.ibb.(λ2 , T1)
	I2 = @. (1.0 - ϵ) * Planck.ibb.(λ2 , T2)
	I_measured =@. two_planck.(λ2 , ϵ , T1 , T2)
end;

# ╔═╡ a7ff8a2d-a81d-4474-b27f-565de2cf5dd3
md""" Cristiansen wavelength: ϵ =$(ϵ[argmax(ϵ)]) at $(λ_max = λ2[argmax(ϵ)]) μm"""

# ╔═╡ 459f54a1-bbf0-4268-8bec-8142d436976a
begin 
	common_kwargs = (; fillrange=0, fillalpha=0.3,dpi=600,legend_background_color=:white, legend_foreground_color = :black, legend_position=:right , grid = true, gridlinewidth=3, gridstyle = :dot,minorgrid=true, box = :on, linewidth = 3)
	
	(xmin , xmax) = extrema(λ_show)
	fl =@. (xmin <= λ2) & (λ2 <= xmax)
	y_lims = extrema(I_measured[fl])
	
	
	em_plot = Plots.plot(λ2 , eint.(λ2)  , label = nothing , title = emissivity_type; common_kwargs...)

	
	
	ppp = Plots.plot(λ2 , I_measured  , label = "sum" ; common_kwargs...)
	Plots.plot!(ppp, λ2 , I1  , label = "Tsample"  ; common_kwargs...)
	Plots.plot!(ppp, λ2 , I2  , label = "Tlamp" ; common_kwargs...)
	
	
	ylabel!(ppp , " I , , W/m²⋅sr⋅μm")
	ylabel!(em_plot , " ϵ")
	ylims!(ppp , y_lims)
	ylims!(em_plot , (0.0,1.0))
	for p in (ppp , em_plot)
		xlabel!(p , " λ , μm ")
		xlims!(p , (xmin , xmax))
		
	end
end

# ╔═╡ 5dafdc88-bd40-4347-aa62-e841d15c1bd7
em_plot

# ╔═╡ 18daa932-fd3a-4056-aa07-4dcf26c7d57a
ppp

# ╔═╡ 0404bf20-57a4-4c7c-bf23-70d3541a6787
begin 

	Tmeas = Dict{Symbol , NamedTuple{(:Tmeas , :ΔT , :ΔTrel) ,Tuple{Float64 , Float64 , Float64}}}()
	for p in pyrometers_vector2
		
	
		is_two_wavelength_pyrometer = Pyrometers.is_spectral_band(p)
		# 
		l_cur = is_two_wavelength_pyrometer ? p.λ : [p.λ[]-0.2,p.λ[]+0.2 ]
		λ_pyr_interp = is_two_wavelength_pyrometer ? collect(range(l_cur...,length=30)) : p.λ

		e_evaluated = eint.(λ_pyr_interp)
		e_averaged = sum(e_evaluated)/length(e_evaluated)

		# setting averaged value of emissivity to pyrometer
		Pyrometers.set_emissivity!(p, e_averaged)
		
		# calculating the measured by the pyrometer value 
		I_cur = two_planck.(λ_pyr_interp , e_evaluated , T1 , T2)
		
		measure_intensity =  is_two_wavelength_pyrometer ? Pyrometers._simpson(λ_pyr_interp , I_cur) : I_cur[]
		# setting emissivity to one
		_tmeasured = Pyrometers.measure(p , measure_intensity, T_starting = T1)
		push!(Tmeas, p.type =>(_tmeasured , _tmeasured - T1 , 100*(_tmeasured - T1)/T1))
	# temperature measured by the current pyrometer 
	end
	
	pyrometers_vector2
end

# ╔═╡ cd9d9742-e9e6-47b9-afae-09e3018e7ebf
Tmeas

# ╔═╡ ac2343b5-6ea5-47b1-9c39-643cdf6d93af
begin 
	kvect = [k for k in keys(Tmeas)]
	(dT, ind_min) = findmin(k->Tmeas[k].ΔT , kvect)
	pyr_smaller_type = kvect[ind_min]
end;

# ╔═╡ 48b184a0-b1df-461e-a3f5-3c8be72ab875
md" Select pyrometer type: $(@bind selected_type Select(kvect , default = pyr_smaller_type) )"

# ╔═╡ 98cfdb54-7338-401f-8dcb-fd7752a70e0e
if use_custom
	p_selected =  Pyrometers.Pyrometer([custom_waves...], type = :C  , ϵ = 1.0)
else
	p_selected = Pyrometers.Pyrometer(selected_type)
end

# ╔═╡ 1ecac1b1-8efb-47d6-9e0a-093487e0c880
p_selected(fix_temperature(PlanckEmitter() , 1234.6))

# ╔═╡ baeafd9c-19bc-4dda-ade2-89b8cf534f88
 begin 
     isurface = p_selected.ϵ[] * Planck.band_power( 1573,  λₗ = custom_waves.left , λᵣ = custom_waves.right) # real temperature is 1700.11
      reflected = (1 - p_selected.ϵ[]) * Planck.band_power( 1873,  λₗ = custom_waves.left , λᵣ = custom_waves.right)
      i = isurface + reflected
 end

# ╔═╡ 1c02bee4-29af-4a5f-8b83-0ee6d5e226be
begin 
	struct MixedRadiationEmitter{E , R , PL1 , PL2}
		e::E
		r::R
		pl1::PL1
		pl2::PL2
		function MixedRadiationEmitter(e)
			r = Pyrometers.SpectralReflectivity(e)
			pl1 = e * Pyrometers.PlanckEmitter()
			pl2 = r * Pyrometers.PlanckEmitter()
			new{typeof(e) , typeof(r) , typeof(pl1) , typeof(pl2)}(e , r , pl1 , pl2)
		end
	end
	(em::MixedRadiationEmitter)(l , t1 , t2) = em.pl1(l , t1) + em.pl2.(l , t2)
	function fix_second_temperature(me::MixedRadiationEmitter , t)
		me.pl1 + Pyrometers.fix_temperature(me.pl2 , t)
	end
end

# ╔═╡ abf9e80e-49e9-4a35-bc68-bb6c4260fa94
begin 
	em1 = MixedRadiationEmitter(e_int_iso)
	em1(2.3 , 1300 , 1800)
	a = fix_second_temperature(em1 , 1800)
	T_measured = p_selected( a , 1200.0 ,  e_int_iso ) # 1200 is the true temperature 
	incident_radiation = Pyrometers.PlanckEmitter()
	Pyrometers.external_source_corrected_temperature(p_selected , T_measured , e_int_iso , 1800.0 , 1.0 , Pyrometers.EnclosureGeometry())
end

# ╔═╡ 5e64a74b-fe14-4461-871b-6b609b5e83cd
plot(λ_show , em1.(λ_show , 300 , 2800))

# ╔═╡ ba559d2d-3bd0-4cd1-8836-8ab0b860c3a5
Pyrometers.integrate(p_selected.λ[1] , p_selected.λ[2]  ,1200.0 ,  a )

# ╔═╡ 4f6c8a96-2347-496a-8d71-1d410fa30ac9
function evaluate_temperature_error(p , e)
	
	T_surf_scan = range(273,1773 , 25)
	T_heater_scan = range(1673, 2973 , 25)

	em_stray_rad = MixedRadiationEmitter(e)

	
	_N = length(T_surf_scan)
	_M = length(T_heater_scan)
	
	ΔT_mat = zeros((_M , _N))
	Tmeas_mat = zeros((_M , _N))
	exclution_error = zeros((_M , _N))
	
	 	for i in 1 : _M
		 t2 = T_heater_scan[i]
		 _a = fix_second_temperature(em_stray_rad , t2) # fixing source temperature returns SpectralQuantity
		Threads.@threads for j in 1 : _N
			t1 =  T_surf_scan[j]
			t_meas =  p( _a , t1 ,  e )
			ΔT_mat[i , j] = (t_meas - t1)/t1
			Tmeas_mat[i , j] = t_meas
			t_stray_excluded = Pyrometers.stray_radiation_corrected_temperature(p , t_meas , e , t2 , 1.0)
			exclution_error[i , j] = (t_stray_excluded - t1)/t1
		end
	end
	return (T_surf_scan, T_heater_scan  , ΔT_mat ,Tmeas_mat ,  exclution_error)
end




# ╔═╡ 10298d52-d411-475f-b7f6-8562ed2a25bc
if is_recalculate 

	(T_surf_scan, T_heater_scan , ΔT_mat ,Tmeas_mat ,  exclution_error)  = evaluate_temperature_error(p_selected , e_int_iso)
	
end

# ╔═╡ 0cb50b02-ab41-416c-8610-c3ff318b117b
if is_recalculate 
	p_res = Plots.surface(T_surf_scan, T_heater_scan , 100.0*ΔT_mat )
end

# ╔═╡ 3e19d251-91f6-4383-bfc8-ffe816570f42
if is_recalculate
	md"""
	Tsurface = $(@bind T_surf_selected Select(T_surf_scan))
	
	Tlamp = $(@bind T_lamp_selected Select(T_heater_scan))
	
	"""
end	

# ╔═╡ ce4f2fdd-16b1-46e8-88a4-a952896b6df8
if is_recalculate

	j_selected = findfirst(t->t == T_surf_selected , T_surf_scan)
	i_selected = findfirst(t->t == T_lamp_selected , T_heater_scan)
	md" Tmeas ± ΔT = $( Tmeas_mat[i_selected , j_selected] ) ± $(ΔT_mat[i_selected , j_selected] *  T_surf_scan[j_selected] )"
end

# ╔═╡ dd1561e2-233f-425a-832f-130b49f0bf0b
if is_recalculate
	out_table = vcat(hcat([0.0] , transpose(T_surf_scan)) , hcat(T_heater_scan, 100*ΔT_mat))
	pretty_table(HTML , out_table)
end

# ╔═╡ 05ec0ea2-cfec-4e91-a46a-68bbdaefd562
if is_recalculate
	out_table_T = vcat(hcat([0.0] , transpose(T_surf_scan)) , hcat(T_heater_scan, Tmeas_mat))
	pretty_table(HTML , out_table_T )
end

# ╔═╡ 805cd62e-a188-4e3f-a870-990530cbc7db
SP = Pyrometers.ScaledPolynomials

# ╔═╡ 5a56e1db-5909-443f-8f8f-e8a640b9bd3e
md"""
	### VI. Multiwavelength pyrometry
	_______________________

	Unlike classical partial radiation pyrometry, which requires setting a constant emissivity in some relatively narrow wavelength region, the **multiwavelength pyrometry**  in theory, allows one to obtain the temperature of a surface without knowing the emissivity. More about multiwavelength pyrometry can be found e.g. in [Multi-spectral pyrometry—a review](https://iopscience.iop.org/article/10.1088/1361-6501/aa7b4b). 

	In order to achieve this goal, the **multiwavelength pyrometry** assumes that the dependence of emissivity on wavelength in some spectral range can be described by some relatively small number (``N``) of parameters. Hence, if you have measured thermal radiation intensity at relatively large number (``M``) of wavelengths, and if ``M>N+1``, you can formulate the optimization problem in space of ``N+1`` optimization variables viz ``\vec{x}= \begin{bmatrix} a_0 , \dots,  a_{N-1} , T\end{bmatrix}``, here `` \begin{bmatrix} a_0 , \dots,  a_{N-1}  \end{bmatrix}`` are ``N``  emissivity approximation coefficients, and the ``(N+1)``'th optimization variable ``T`` is the temperature.	The **multiwavelength pyrometry** optimization problem has several features that can be utilized in order to obtain a computationally-effective algorithm:

	* First, the emissivity approximation is a linear problem, this means that emissivity approximation coefficients can be taken independent of both the optimization variables and the independent variables (wavelength)
	* Second, a highly non-linear term (viz Planck function) depends on only one of the optimizaiton variables
	* And third, the target function is the product of linear and non-linear optimization problems

	Mathematical consequencies of these features are described in this repository supplementary materials [download pdf](https://manarom.github.io/BandPyrometry.jl/assets/supplementary_v_0_0_1.pdf) in more details.

	In `Pyrometers.jl` , the real surface thermal emission is approximated as a product of Planck function (ideal surface thermal emission) and linear (with respect to the optimization variables e.g. polynomial coefficients) approximation of spectral emissivity.	

	In `Pyrometers.jl` the spectral emissivity is approximated as a linear combination of basis functions:
	"""

# ╔═╡ 253847d7-87fa-461d-8361-fac6c6facf73
md"""
``\epsilon(λ)=\sum_{n=0}^{N-1}a_n \cdot \phi_n(λ)``

where ``\phi_n`` is the basis function column vector, e.g. for standard basis it is:

``\mathbf{ \phi_n(λ)  = \begin{bmatrix} λ₁ⁿ  \\ \vdots \\ λ_{M}ⁿ \end{bmatrix}
}``

Full emission spectrum of a real surface is calculated as:

``\mathbf{ 
I_{real\:surface} (\lambda,T) =I_{blackbody}(\lambda,T) \cdot  \epsilon(λ)=I_{blackbody}(\lambda,T) \cdot \sum_{n=0}^{N-1}a_n \cdot \phi_n(λ)
} ``

Now, the optimization problem can be formulated:

``\vec{x}^*=argmin\{F(\vec{x})\}`` 

``F(\vec{x})=\sum_{i=1}^{M}[y_i - I_{blackbody}(\lambda_i,T) \cdot (\sum_{n=0}^{N-1}a_n \cdot \phi_n(λ_i))]^2`` 

``\vec{x} = [\vec{a},T]^t``

where ``\vec{a}`` is the column vector of emissivity approximation (  ``[]^t`` stays for transposition), ``\vec{x}^*`` is the local minimum

To solve this optimization problem `Pyrometers.jl` package provides functions to evaluate the discrepancy function ``F``, ``\nabla F`` and ``\nabla^2F`` which are needed to solve the optimization problem using zero, first or second order optimization algorithms. It also has special type to work with the emissivity linear approximation.

"""

# ╔═╡ 86b4b811-7c70-49a5-92f2-4ba409a0ef32
md"""

####  Emissivity approximation functions
_______________________

To approximate the emissivity `Pyrometers.jl`  uses several polynomial bases:
(from [`ScaledPolynomials.jl`](https://https://github.com/Manarom/ScaledPolynomials.jl) package)
* Standard  
* Chebyshev  
* Legendre 
* Trigonometric basis: ``\phi_n(λ) = sin(\pi \cdot n \lambda)``  for odd n, and  ``\phi_n(λ) = cos(\pi \cdot n \lambda)`` for even n
* Bernstein basis of degree ``D``: ``\phi_{k}^{n}(\lambda) = (\begin{matrix} n \\ k \end{matrix}) (\frac {\lambda - a}{b - a})^k(\frac {b-\lambda}{b-a})^{n-k}`` - Bernstein basis function for ``\lambda \in [\lambda_a...\lambda_b]`` , ``k \in [0...n]`` , ``a = -1``, ``b = 1``

All basis vectors are stored in a structure called `VanderMatrix`, this type: 

1. Stores basis vectors for selected polynomial type and degreee  - columns of matrix ``V``: ``V= \begin{bmatrix} \vec{\phi_1} , \dots,  \vec{\phi_n} \end{bmatrix}``

2. The resulting emissivity can be calculated as a product of `VanderMatrix` and the vector of emissivity approximation polynomial coefficients vector: ``\vec{\epsilon}=V\cdot\vec{a}``

"""

# ╔═╡ 4894c2ab-4db5-4b7c-9c4a-9dd1c5345c28
@bind  λ_fit_vand PlutoUI.combine() do Child
	md"""
	Emissivity fitting spectral range, ``\mu m`` : \
	``\lambda_{left}`` = $(
		Child(Slider(0.1:0.1:16,default=4.0,show_value = true))
	)   -- 
	 $(
		Child(Slider(0.1:0.1:16,default=8.0,show_value = true))
	)  ``\lambda_{right}`` 
	"""
end

# ╔═╡ c2582621-54fb-44b1-a3ea-4cfacb6062ff
md"Select polynomial basis type $(@bind em_approx_poly_type Select(collect(keys(Pyrometers.ScaledPolynomials.SUPPORTED_POLYNOMIAL_TYPES)),default = :bernsteinsym))"

# ╔═╡ c7554489-1d97-4b9a-a3e9-4be84c82b552
md"Set polynomial degree : $(@bind poly_fit_degree Select(0:10,default=3)) (the polynomial degree = numer of basis functions - 1, thus, zero order polynomial is an all-units vector)"

# ╔═╡ 5815a317-233a-493a-a8be-03dc7d608c0e
begin 
	N1 = 25
	λ_fit_vec = collect(range(λ_fit_vand...,length=N1)) 
	
	interpolated_values =rt_emissivity_interpolation(λ_fit_vec) # interpolating data at λ points
	PolyType = Pyrometers.ScaledPolynomials.SUPPORTED_POLYNOMIAL_TYPES[em_approx_poly_type]{poly_fit_degree + 1 , Float64}
	Vander_test = Pyrometers.ScaledPolynomials.VanderMatrix(SVector{N1}(λ_fit_vec),PolyType()) # creating new matrix 
	(a_fit_check,fitted_value,goodness_of_fit) = Pyrometers.ScaledPolynomials.polyfit(Vander_test,λ_fit_vec,interpolated_values)# fitting polynomial coefficients
	plot(rt_emissivity_data[:,1],rt_emissivity_data[:,2],label = "ϵ real")
	scatter!(λ_fit_vec,fitted_value, label="ϵ fitted"; plot_common_args...)
	xlabel!("Wavelength, μm")
	ylabel!("Emissivity")
end

# ╔═╡ 77b352db-6734-4501-b9d3-f63ff2adbaf7
pretty_table(HTML,hcat(["a$(i)" for i in 0:1:poly_fit_degree],a_fit_check), top_left_string="Table of the coefficients of emissivity linear approximation in the band from $(λ_fit_vand[1]) to $(λ_fit_vand[2]) μm using $(em_approx_poly_type) bases type,  the goodness of fit = $(goodness_of_fit)"  , column_labels = ["coeff","val"])

# ╔═╡ bf58aa17-dc84-4bee-94d3-25895b12f553
md"Measured spectrum fitting region:"

# ╔═╡ c8647683-a25e-4c04-bae0-52a5f40233e9
md"Select the polynomial type = $(@bind poly_type confirm(Select(collect(keys(Pyrometers.ScaledPolynomials.SUPPORTED_POLYNOMIAL_TYPES)),default = :bernsteinsym)))"

# ╔═╡ 7957b928-29db-4342-9993-15b63023883b
md"Set polynomial degree for real emissivity: $(@bind real_poly_degree confirm(Select(0:6,default=3))) (the polynomial degree= numer of basis functions-1, thus zero order polynomial is constant)"

# ╔═╡ 3b01166a-451e-4e15-ae34-7049703331f2
md"""
The following two figures show:

1)basis vectors for selected polynomial (columns of `VanderMatrix.v`): ``V``

2)the resulting emissivity, calculated as a product of `VanderMatrix` and the vector of emissivity approximation polynomial: ``\vec{\epsilon}= \begin{bmatrix} \vec{\phi_1} , \dots,  \vec{\phi_n} \end{bmatrix}\cdot\vec{a} = V\cdot\vec{a}``

"""

# ╔═╡ 6a386feb-48a9-40ac-8fac-be492183ed2b
@bind  a_real PlutoUI.combine() do Child
	md"""
	a0 = $(
		Child(Slider(-1:0.01:1,default=0.6,show_value = true))
	) \
	a1 = $(
		Child(Slider(-1:0.01:1,default=0.8,show_value = true))
	)\
	a2 = $(
		Child(Slider(-1:0.01:1,default=0.8,show_value = true))
	)\
	a3 = $(
		Child(Slider(-1:0.01:1,default=0.8,show_value = true))
	)\
	a4 = $(
		Child(Slider(-1:0.01:1,default=0.8,show_value = true))
	)\
	a5 = $(
		Child(Slider(-1:0.01:1,default=0.8,show_value = true))
	)\
	"""
end

# ╔═╡ 10c7b1f0-3565-456d-a4bf-84aff0aca60b
begin 
	poly_obj = Pyrometers.ScaledPolynomials.SUPPORTED_POLYNOMIAL_TYPES[poly_type](a_real[1:real_poly_degree])
	scaled_poly = SP.ScaledPolynomial(poly_obj , xmin = λ_fit_vand[1] , xmax =λ_fit_vand[2])
	λ = MVector{50}(range(λ_fit_vand[1] , λ_fit_vand[2] , 50))
end

# ╔═╡ 62a86ff0-98a9-4acf-bed0-91a5eab24209
md" Ttrue = $(@bind Ttrue Slider(100:1.0:3000 , default = 1000.0 , show_value = true))"

# ╔═╡ 646705d8-42e9-4204-ac7d-f2d63425b63c
md" ϵ bounds = $(@bind e_bounds RangeSlider(0.01:1e-2:1.0))"

# ╔═╡ 3bc1a1a0-5d6a-405a-ac3a-e67541e023c5
md" Tstarting $(@bind Tstarting Slider(100:1.0:3000 , default = 1200.0 , show_value = true))"

# ╔═╡ c0f25834-2bc2-4c65-aefa-466d8017d461
begin 
	p_em = plot(λ, scaled_poly.(λ),label=nothing,linewidth=6; plot_common_args...)
	title!(raw"""Generated "real" surface spectral emissivity""")
    xlabel!("Wavelength, μm")
	ylabel!("ϵ")
	p_em
end

# ╔═╡ 229e5f0b-5d91-4996-9ca8-beaa689ea2da
@bind refit Button("refit")

# ╔═╡ ac1babcb-5dad-4bc6-bb4c-c707dbd57fa0
ParticleSwarm

# ╔═╡ 963fb46a-0ea3-48b5-b62f-e37c8fde1864
begin 
	e_surf = IsothermalSpectralQuantity(scaled_poly)
	i_measured = Pyrometers.fix_temperature(e_surf * PlanckEmitter() , Ttrue)
	multiwavelength_pyro = Pyrometers.MultiWavelengthPyrometer{length(λ) , poly_fit_degree , :bernstein}(λ , i_measured=i_measured)
end

# ╔═╡ f6c0cbdf-d2a6-47c7-bc58-edc071760df9
begin 
	refit 
	multiwavelength_pyro(emissivity_range = extrema(e_bounds) , optimizer = GradientDescent)

end

# ╔═╡ a5f30290-343d-4a3a-ad85-589fbe8ed570
begin 
	refit
	e_fitted = emissivity_poly(multiwavelength_pyro)
	plot(λ , e_surf.(λ) , label = "true : T=$(Ttrue)")
	plot!(λ , e_fitted.(λ) , label = "fitted, T=$(Pyrometers.temperature(multiwavelength_pyro))")
	title!("Emissivity identification result and measured temperature")
end

# ╔═╡ 6c712eb2-8e41-4d36-a4e3-077905eb4214
begin 
	sence = Pyrometers.sensitivity(multiwavelength_pyro)
	sence_plot = plot(;plot_common_args...)
	for (i , c) in enumerate(eachcol(sence.S))
		plot!(sence.l , c , label = i)
	end
	sence_plot
end

# ╔═╡ Cell order:
# ╟─30743a02-c643-4bdc-837e-b97299f9520a
# ╠═5e712312-0fc7-4205-84cc-834d57b814a3
# ╠═f728a59d-c78c-45af-a9e5-656be490eb4f
# ╟─abdc809b-b53c-4dff-ba6f-c636c73f3fca
# ╠═0bee71f4-5961-4c80-8592-2b6c0d1b58a8
# ╠═bf833e74-f9e7-4b60-b6bc-2a6a58c5c901
# ╟─05c05c84-02d4-4b7f-83df-bd1fa3e4ee4d
# ╟─171409eb-22b5-4bc5-a8e2-eac0932a24f3
# ╟─643d9ff3-3a09-46c9-9013-92d111ccb229
# ╟─d5ee3913-66be-47d7-a755-699ba64b4f98
# ╟─d442014a-20e6-4be4-ac7f-f13de329dec5
# ╟─27b3c586-9eb0-4a51-b9ca-a9c0379fccdf
# ╟─f22d22b6-5d98-4cc4-998f-a53e92809618
# ╟─b9bee300-59a4-4c7a-b525-439f5c62253e
# ╟─7cc110e9-7655-4dfc-b1e0-ab3905866425
# ╟─4d6337aa-cfc7-4154-a395-5aa53e23d01a
# ╟─03d76e64-ebf4-432b-b9be-d4cb26275f55
# ╟─59b01938-9235-4d39-be4e-e4fdf9f94418
# ╟─6e3f01fb-a0bb-438b-a95a-2c503cc4961e
# ╟─b7117f7b-ad22-4b09-8aa3-3f6a3878c149
# ╟─5efd8d48-7fff-48fd-b361-2f06cd8bad53
# ╟─af619a4d-b700-4996-960c-c9d0e75eac6a
# ╟─0aad4eb6-e5dc-40e3-82b6-4d64f5b1af4f
# ╟─e9e16b2b-88d7-40b6-997f-3a40e673faa8
# ╟─2ad3ec82-54a2-49ac-94ef-579f808dfb1a
# ╟─63d5b5bb-5e7d-4e83-83ae-3972f2465e28
# ╟─95baff6c-2af3-401d-9db3-d6e6e31ab70a
# ╟─27203e6b-9283-4da3-beb2-9f671100860f
# ╟─5010ec7e-b83c-47c6-8c84-17bbcb8e40a1
# ╟─8a066ee5-80e9-462f-9a61-15851468aa63
# ╟─b7fac177-c211-4635-992f-e6473be7bdae
# ╟─e2a9aa39-2490-4681-89d3-a01f058f6feb
# ╟─02eee968-ff43-4b82-8d68-efede1a220dd
# ╠═1ecac1b1-8efb-47d6-9e0a-093487e0c880
# ╟─e9d0e216-07b5-4e41-925a-0c530e6ddf7b
# ╟─9b08b767-7e8f-4483-9f2f-226022ce10e4
# ╟─144b40ea-71c7-421f-8117-eab267ea5daf
# ╟─c69acbf6-94fb-4ac3-8d56-d1f9dda11440
# ╟─a861d56f-f6c9-4754-b7a9-ed63713f1f2f
# ╟─7071a6f4-e296-4e53-8e6f-24f1f038c1a5
# ╟─712828a7-fb54-42e6-95fc-233243190f59
# ╟─f763d449-2a7a-4008-a183-823a774bc25e
# ╟─667f7c30-56e0-461f-b35b-c924007eb9f2
# ╟─d08ec8f2-7689-4043-9f28-da06ab0124b9
# ╟─c5ac80ee-8143-4c28-bbff-2ac761c71fac
# ╟─0c9fe7b1-374c-4fb8-9cfe-9337389713bf
# ╟─bc2d93ae-6c30-462c-96e0-30fdb84d7c63
# ╟─d3199b6e-9779-4def-b701-fe85d1035045
# ╟─8a6fe87d-0f7c-4577-9ac2-ea1ecf71016b
# ╟─72947d97-0a97-4064-a2fe-08d19dec0f0e
# ╟─9ce196c1-8915-46da-9aba-f13d7655959a
# ╟─3581aa29-714b-422a-8feb-d1a0c3ebeec7
# ╟─86af6afc-b28a-4e84-952a-bd29710374f8
# ╟─efc35420-d0e6-4795-94b6-d43289b4de44
# ╟─6342e92b-4434-4e4b-aa2f-56405277caed
# ╟─f181980f-bf72-4468-8daa-9461c6c901e0
# ╟─5dafdc88-bd40-4347-aa62-e841d15c1bd7
# ╟─a7ff8a2d-a81d-4474-b27f-565de2cf5dd3
# ╠═18daa932-fd3a-4056-aa07-4dcf26c7d57a
# ╟─36ba2396-bb5e-4d22-a58e-9ab27cd18b2d
# ╟─91bbd553-4e4a-431d-9d54-b0f4882fd426
# ╟─15f1519b-d924-4fa9-b212-eebba75c544a
# ╟─cd9d9742-e9e6-47b9-afae-09e3018e7ebf
# ╟─0404bf20-57a4-4c7c-bf23-70d3541a6787
# ╟─8cef05a1-2974-4c38-b73e-fa706f347fcc
# ╟─10298d52-d411-475f-b7f6-8562ed2a25bc
# ╟─459f54a1-bbf0-4268-8bec-8142d436976a
# ╠═1811e43c-f7db-47b1-9b83-bb38455d7db3
# ╟─ac2343b5-6ea5-47b1-9c39-643cdf6d93af
# ╟─54339700-71fd-48bf-a2ef-0c3267b9d81b
# ╟─7f76bc22-77c3-4eb3-9d49-588e653df2e7
# ╟─a0197a9a-34bf-4a3e-af8a-c23ea777f482
# ╟─7793976e-6714-4bc1-9d97-412dd2a67480
# ╟─48b184a0-b1df-461e-a3f5-3c8be72ab875
# ╟─98cfdb54-7338-401f-8dcb-fd7752a70e0e
# ╟─baeafd9c-19bc-4dda-ade2-89b8cf534f88
# ╟─0cb50b02-ab41-416c-8610-c3ff318b117b
# ╟─3e19d251-91f6-4383-bfc8-ffe816570f42
# ╟─ce4f2fdd-16b1-46e8-88a4-a952896b6df8
# ╟─dd1561e2-233f-425a-832f-130b49f0bf0b
# ╟─05ec0ea2-cfec-4e91-a46a-68bbdaefd562
# ╠═e480137d-b6d9-4e18-92f0-640292bbb5f0
# ╠═1c02bee4-29af-4a5f-8b83-0ee6d5e226be
# ╠═abf9e80e-49e9-4a35-bc68-bb6c4260fa94
# ╟─5e64a74b-fe14-4461-871b-6b609b5e83cd
# ╟─ba559d2d-3bd0-4cd1-8836-8ab0b860c3a5
# ╟─4f6c8a96-2347-496a-8d71-1d410fa30ac9
# ╠═805cd62e-a188-4e3f-a870-990530cbc7db
# ╟─5a56e1db-5909-443f-8f8f-e8a640b9bd3e
# ╟─253847d7-87fa-461d-8361-fac6c6facf73
# ╟─86b4b811-7c70-49a5-92f2-4ba409a0ef32
# ╟─4894c2ab-4db5-4b7c-9c4a-9dd1c5345c28
# ╟─c2582621-54fb-44b1-a3ea-4cfacb6062ff
# ╟─c7554489-1d97-4b9a-a3e9-4be84c82b552
# ╟─5815a317-233a-493a-a8be-03dc7d608c0e
# ╟─77b352db-6734-4501-b9d3-f63ff2adbaf7
# ╟─bf58aa17-dc84-4bee-94d3-25895b12f553
# ╟─c8647683-a25e-4c04-bae0-52a5f40233e9
# ╟─7957b928-29db-4342-9993-15b63023883b
# ╟─3b01166a-451e-4e15-ae34-7049703331f2
# ╟─10c7b1f0-3565-456d-a4bf-84aff0aca60b
# ╟─6a386feb-48a9-40ac-8fac-be492183ed2b
# ╟─62a86ff0-98a9-4acf-bed0-91a5eab24209
# ╠═646705d8-42e9-4204-ac7d-f2d63425b63c
# ╟─3bc1a1a0-5d6a-405a-ac3a-e67541e023c5
# ╟─c0f25834-2bc2-4c65-aefa-466d8017d461
# ╠═229e5f0b-5d91-4996-9ca8-beaa689ea2da
# ╠═f6c0cbdf-d2a6-47c7-bc58-edc071760df9
# ╠═ac1babcb-5dad-4bc6-bb4c-c707dbd57fa0
# ╠═a5f30290-343d-4a3a-ad85-589fbe8ed570
# ╠═6c712eb2-8e41-4d36-a4e3-077905eb4214
# ╠═963fb46a-0ea3-48b5-b62f-e37c8fde1864
# ╠═a8a27a41-95a2-41eb-a1cf-d3d51b2ec52e
