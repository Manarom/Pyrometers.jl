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
	using QuadGK
	src_dir = joinpath(abspath(joinpath(notebook_dir,"..")),"src")
end;

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

# ╔═╡ abdc809b-b53c-4dff-ba6f-c636c73f3fca
const Planck = Pyrometers.Planck

# ╔═╡ ba23c985-74c4-41f3-8bc4-f7287e30e47f
#using Revise,StaticArrays,OrderedCollections,Optimization,OptimizationOptimJL,LaTeXStrings,Interpolations,Plots,PlutoUI,DelimitedFiles , ForwardDiff , Roots , QuadGK

# ╔═╡ 9cd8fe6d-dcf9-472e-a019-19b4c1a182ed
#includet(joinpath(src_dir,"RadiationPyrometers.jl"))

# ╔═╡ 6535717c-99ae-4e8e-94aa-d600f02537ed
# ╠═╡ disabled = true
#=╠═╡
 Planck = RadiationPyrometers.Planck
  ╠═╡ =#

# ╔═╡ 171409eb-22b5-4bc5-a8e2-eac0932a24f3
PlutoUI.TableOfContents(indent=true, depth=4, aside=true)

# ╔═╡ d5ee3913-66be-47d7-a755-699ba64b4f98
md"""
### Introduction

This notebook demonstrates two examples of using the **RadiationPyrometers.jl** package. The purpose of this small package is to create a virtual pyrometer that can be used to calculate the emissivity of a real-life pyrometer, enabling the measured temperature to be adjusted to match the actual temperature of the heated object.  

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

# ╔═╡ 8a066ee5-80e9-462f-9a61-15851468aa63
md"""
### II. Partial radiation pyrometry
_______________________

As far as the `blackbody` thermal radiation energy strongly depends on temperature, this quantity can be used to measure the temperature of a real surface. This is the general idea of partial radiation pyrometry: **measure intensity to get the temperature**. As far as the intensity is a directional quantity, a pyrometer needs collimating optics (a telescope).The real surfaces emissivity often varies sufficiently with the wavelength, at the same time, partial radiation pyrometers assume constant emissivity (so-called `grey`-band approximation). Thus, for industrial purposes, it is useful to have several pyrometers, each working within a relatively narrow spectral band. In the spectral range of a partial radiation pyrometer, emissivity should not vary significantly to make the assumption of constant emissivity relevant.

The  [RadiationPyrometers.jl](https://manarom.github.io/RadiationPyrometers.jl) package provides several function to work with `virtual` partial radiation pyrometers.
"""

# ╔═╡ b7fac177-c211-4635-992f-e6473be7bdae
md"""
Dictionary **RadiationPyrometers.DefaultPyrometersTypes** contains default pyrometers type names together with spectral range. Custom pyrometer can be created by providing its type (name), wavelength or wavelentgh range, working emissivity: 
``` julia
# creating custom pyrometer objects
p_custom = RadiationPyrometers.Pyrometer(type = "Custom",λ = [2.5, 3.7],ϵ=0.65)
```
"""

# ╔═╡ e2a9aa39-2490-4681-89d3-a01f058f6feb
pretty_table(HTML,Pyrometers.DefaultPyrometersTypes,top_left_string ="Table of default pyrometers types provied by `RadiationPyrometers.jl` package and corresponding wavelength regions",wrap_table_in_div=true)


# ╔═╡ 02eee968-ff43-4b82-8d68-efede1a220dd
md"""
After creating the **Pyrometer** object, it can be used to "measure" the temperature viz convert the intensity of a real surface to its temperature, according to the pyrometer's spectral range using **RadiationPyrometers.measure(pyr,measured_intensity)** function (the intensity should be provided in correct units [W/m²⋅sr⋅μm]). It is quite in practice to know the real temperature of the surface at some point, in this case the emissivity of virtual pyrometer can be adjusted to make the "measured" temperature be equal to the real one. This can be done using 
```julia
	RadiationPyrometers.fit_ϵ!(p::Pyrometer,Tmeasured::Float64,Treal::Float64) 
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
### III. Blackbody reference source emissivity
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
 ### IV Two temperature emission
"""

# ╔═╡ 9ce196c1-8915-46da-9aba-f13d7655959a
md"""
 	Now the emission is proportional to the spectral intensity of two mixed BB sources with different temperatures.

``I_{meas} = \epsilon\mathcal{b}(\lambda , T_1) + (1 - \epsilon)\mathcal{b}(\lambda , T_2)``

Check is it is possible to measure the surface temperature in a presence of external emission with much higher temperature
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
		ϵ = ϵ_fixed
		eint = Returns(ϵ_fixed)
	else
		_data = readdlm(joinpath(data_folder , emissivity_type))
		if contains(emissivity_type , "quarz")
			@. _data[:,2] = 1.0 - _data[:,2]
		end
		eint = linear_interpolation(_data[:,1] , _data[:,2] , extrapolation_bc=Line())
		ϵ = eint.(λ2)
		 
	end
end;

# ╔═╡ a7ff8a2d-a81d-4474-b27f-565de2cf5dd3
md""" Cristiansen wavelength: ϵ =$(ϵ[argmax(ϵ)]) at $(λ_max = λ2[argmax(ϵ)]) μm"""

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
md" ### Recalculate T matrix $(@bind is_recalculate CheckBox(false))"

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

λleft = $(Child("left", NumberField(0.1:1e-3:20 , default = 7.2)))

λright = $(Child("right", NumberField(0.1:1e-3:20 , default = 7.3)))

"""
end

# ╔═╡ 969907ad-0c38-4dfb-8e2d-ac25619fbe5f
custom_waves

# ╔═╡ 5b647454-e5fc-4c65-8a0a-8eb499955cc1
if is_recalculate && use_custom
	p_custom = Pyrometers.Pyrometer([custom_waves...],type = :C  , ϵ = 1.0)
end

# ╔═╡ c9634225-fa52-4747-8f45-3511141bd164
md" #### Try plotly! : $(@bind is_use_plotly CheckBox(false))"

# ╔═╡ e480137d-b6d9-4e18-92f0-640292bbb5f0
function two_planck(l , ϵ , T1 , T2)
	return ϵ * Planck.ibb(l , T1) + (1.0 - ϵ) * Planck.ibb(l , T2)
end

# ╔═╡ 6342e92b-4434-4e4b-aa2f-56405277caed
begin 
	I1 = @. ϵ * Planck.ibb.(λ2 , T1)
	I2 = @. (1.0 - ϵ) * Planck.ibb.(λ2 , T2)
	I_measured =@. two_planck.(λ2 , ϵ , T1 , T2)
end;

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

# ╔═╡ 10298d52-d411-475f-b7f6-8562ed2a25bc
if is_recalculate 
	
	T1_scan = 273.0:50:1373
	T2_scan = 1273.0:50:3800


	if use_custom
		p_selected = p_custom
	else
		p_selected = filter(p->p.type == selected_type , pyrometers_vector2)[]
	end
	
	_is2 = Pyrometers.is_spectral_band(p_selected)
		# 
	l = _is2 ? p_selected.λ : [p_selected.λ[]-0.2,p_selected.λ[]+0.2 ]
	
	_λ_pyr_interp = _is2 ? collect(range(l...,length=30)) : p_selected.λ

	e_s = eint.(_λ_pyr_interp)
	e_avg = Pyrometers.integral_emissivity(p_selected , eint , 1273.15)
	
	#e_avg = sum(e_s)/length(e_s)

	# setting averaged value of emissivity to pyrometer
	Pyrometers.set_emissivity!(p_selected, e_avg)
	
	# calculating the measured by the pyrometer value 
	
	_N = length(T1_scan)
	_M = length(T2_scan)
	
	ΔT_mat = zeros((_M , _N))
	Tmeas_mat = zeros((_M , _N))
	
	 	for j in 1 : _N
		for i in 1 : _M
			t1 =  T1_scan[j]
			t2 = T2_scan[i]
			i_measured_selected = two_planck.(_λ_pyr_interp , e_s , t1, t2)
			
			i_measured_int =  _is2 ? quadgk(l-> two_planck(l , eint(l) , t1, t2) , p_selected.λ[1] , p_selected.λ[2] )[1] : i_measured_selected[]

			t_meas = try 
					p_selected(i_measured_int)
				catch 
					0.0
			end
			
			ΔT_mat[i , j] = (t_meas - t1)/t1
			Tmeas_mat[i , j] = t_meas
		end
	end
	
end

# ╔═╡ fe1b4c33-1be4-46ab-ba70-b58504e169a6
dfff = quadgk(l-> two_planck(l , eint(l) , 1573 , 1873) , p_selected.λ[1] , p_selected.λ[2])[1]

# ╔═╡ baeafd9c-19bc-4dda-ade2-89b8cf534f88
 begin 
     isurface = p_selected.ϵ[] * Planck.band_power( 1573,  λₗ = custom_waves.left , λᵣ = custom_waves.right) # real temperature is 1700.11
      reflected = (1 - p_selected.ϵ[]) * Planck.band_power( 1873,  λₗ = custom_waves.left , λᵣ = custom_waves.right)
      i = isurface + reflected
 end

# ╔═╡ b82e1650-ff56-45b1-baa1-ebf534400d73
Pyrometers.integral_emissivity(p_selected , eint , 1573.15)

# ╔═╡ cf6322f0-0802-4261-b78f-a0f97a8ae2ad
Pyrometers.set_emissivity!(p_selected , Pyrometers.integral_emissivity(p_selected , eint , 1500.0))

# ╔═╡ 7450bf72-3a21-4322-995b-36e16e62fd18
T_measured = p_selected(i)

# ╔═╡ 43576cb6-87d6-4bcc-90bd-a7ab754e1245
Pyrometers.corrected_temperature(p_selected ,  T_measured , 1873 , 1.0)

# ╔═╡ ca361b35-22b3-4c8f-a1c2-ca5ffc594eff
p_selected

# ╔═╡ 0cb50b02-ab41-416c-8610-c3ff318b117b
if is_recalculate 
	if is_use_plotly
		tr2 = PlutoPlotly.surface(x = T2_scan , y = T1_scan, z=100.0 * ΔT_mat, colorscale="Viridis")
		
		layout2 = Layout(
    		width=800, 
    		height=600, 
    		autosize=true,
    		margin=attr(l=0, r=0, b=0, t=50),  # Minimize margins
    		scene=attr(
        		yaxis_title="T образца",
            	xaxis_title="T лампы", 
            	zaxis_title="100%*ΔT/T1"

    		)
		)
		p_res = PlutoPlotly.plot(tr2 , layout2)
	else
		p_res = Plots.surface(T1_scan, T2_scan , 100.0*ΔT_mat )
	end
end

# ╔═╡ 3e19d251-91f6-4383-bfc8-ffe816570f42
if is_recalculate
	md"""
	Tsurface = $(@bind T_surf_selected Select(T1_scan))
	
	Tlamp = $(@bind T_lamp_selected Select(T2_scan))
	
	"""
end	

# ╔═╡ ce4f2fdd-16b1-46e8-88a4-a952896b6df8
if is_recalculate

	j_selected = findfirst(t->t == T_surf_selected , T1_scan)
	i_selected = findfirst(t->t == T_lamp_selected , T2_scan)
	md" Tmeas ± ΔT = $( Tmeas_mat[i_selected , j_selected] ) ± $(ΔT_mat[i_selected , j_selected] *  T1_scan[j_selected] )"
end

# ╔═╡ dd1561e2-233f-425a-832f-130b49f0bf0b
if is_recalculate
	out_table = vcat(hcat([0.0] , transpose(T1_scan)) , hcat(T2_scan, 100*ΔT_mat))
	pretty_table(HTML , out_table)
end

# ╔═╡ 05ec0ea2-cfec-4e91-a46a-68bbdaefd562
if is_recalculate
	out_table_T = vcat(hcat([0.0] , transpose(T1_scan)) , hcat(T2_scan, Tmeas_mat))
	pretty_table(HTML , out_table_T )
end

# ╔═╡ Cell order:
# ╠═30743a02-c643-4bdc-837e-b97299f9520a
# ╠═5e712312-0fc7-4205-84cc-834d57b814a3
# ╠═abdc809b-b53c-4dff-ba6f-c636c73f3fca
# ╠═ba23c985-74c4-41f3-8bc4-f7287e30e47f
# ╠═9cd8fe6d-dcf9-472e-a019-19b4c1a182ed
# ╟─6535717c-99ae-4e8e-94aa-d600f02537ed
# ╟─171409eb-22b5-4bc5-a8e2-eac0932a24f3
# ╟─d5ee3913-66be-47d7-a755-699ba64b4f98
# ╟─d442014a-20e6-4be4-ac7f-f13de329dec5
# ╟─27b3c586-9eb0-4a51-b9ca-a9c0379fccdf
# ╟─f22d22b6-5d98-4cc4-998f-a53e92809618
# ╟─b9bee300-59a4-4c7a-b525-439f5c62253e
# ╟─7cc110e9-7655-4dfc-b1e0-ab3905866425
# ╟─4d6337aa-cfc7-4154-a395-5aa53e23d01a
# ╟─03d76e64-ebf4-432b-b9be-d4cb26275f55
# ╟─2ad3ec82-54a2-49ac-94ef-579f808dfb1a
# ╟─8a066ee5-80e9-462f-9a61-15851468aa63
# ╟─b7fac177-c211-4635-992f-e6473be7bdae
# ╟─e2a9aa39-2490-4681-89d3-a01f058f6feb
# ╟─02eee968-ff43-4b82-8d68-efede1a220dd
# ╟─9b08b767-7e8f-4483-9f2f-226022ce10e4
# ╟─144b40ea-71c7-421f-8117-eab267ea5daf
# ╟─c69acbf6-94fb-4ac3-8d56-d1f9dda11440
# ╟─a861d56f-f6c9-4754-b7a9-ed63713f1f2f
# ╠═7071a6f4-e296-4e53-8e6f-24f1f038c1a5
# ╟─712828a7-fb54-42e6-95fc-233243190f59
# ╟─f763d449-2a7a-4008-a183-823a774bc25e
# ╟─667f7c30-56e0-461f-b35b-c924007eb9f2
# ╟─d08ec8f2-7689-4043-9f28-da06ab0124b9
# ╟─c5ac80ee-8143-4c28-bbff-2ac761c71fac
# ╟─0c9fe7b1-374c-4fb8-9cfe-9337389713bf
# ╟─bc2d93ae-6c30-462c-96e0-30fdb84d7c63
# ╠═d3199b6e-9779-4def-b701-fe85d1035045
# ╠═8a6fe87d-0f7c-4577-9ac2-ea1ecf71016b
# ╟─72947d97-0a97-4064-a2fe-08d19dec0f0e
# ╟─9ce196c1-8915-46da-9aba-f13d7655959a
# ╟─3581aa29-714b-422a-8feb-d1a0c3ebeec7
# ╟─86af6afc-b28a-4e84-952a-bd29710374f8
# ╟─efc35420-d0e6-4795-94b6-d43289b4de44
# ╟─6342e92b-4434-4e4b-aa2f-56405277caed
# ╟─f181980f-bf72-4468-8daa-9461c6c901e0
# ╟─5dafdc88-bd40-4347-aa62-e841d15c1bd7
# ╟─a7ff8a2d-a81d-4474-b27f-565de2cf5dd3
# ╟─18daa932-fd3a-4056-aa07-4dcf26c7d57a
# ╟─36ba2396-bb5e-4d22-a58e-9ab27cd18b2d
# ╟─91bbd553-4e4a-431d-9d54-b0f4882fd426
# ╟─15f1519b-d924-4fa9-b212-eebba75c544a
# ╟─cd9d9742-e9e6-47b9-afae-09e3018e7ebf
# ╟─0404bf20-57a4-4c7c-bf23-70d3541a6787
# ╟─8cef05a1-2974-4c38-b73e-fa706f347fcc
# ╟─459f54a1-bbf0-4268-8bec-8142d436976a
# ╟─1811e43c-f7db-47b1-9b83-bb38455d7db3
# ╟─ac2343b5-6ea5-47b1-9c39-643cdf6d93af
# ╟─54339700-71fd-48bf-a2ef-0c3267b9d81b
# ╟─7f76bc22-77c3-4eb3-9d49-588e653df2e7
# ╟─a0197a9a-34bf-4a3e-af8a-c23ea777f482
# ╠═7793976e-6714-4bc1-9d97-412dd2a67480
# ╠═969907ad-0c38-4dfb-8e2d-ac25619fbe5f
# ╟─48b184a0-b1df-461e-a3f5-3c8be72ab875
# ╠═5b647454-e5fc-4c65-8a0a-8eb499955cc1
# ╠═fe1b4c33-1be4-46ab-ba70-b58504e169a6
# ╠═baeafd9c-19bc-4dda-ade2-89b8cf534f88
# ╠═b82e1650-ff56-45b1-baa1-ebf534400d73
# ╠═cf6322f0-0802-4261-b78f-a0f97a8ae2ad
# ╠═7450bf72-3a21-4322-995b-36e16e62fd18
# ╠═43576cb6-87d6-4bcc-90bd-a7ab754e1245
# ╠═ca361b35-22b3-4c8f-a1c2-ca5ffc594eff
# ╠═10298d52-d411-475f-b7f6-8562ed2a25bc
# ╟─c9634225-fa52-4747-8f45-3511141bd164
# ╠═0cb50b02-ab41-416c-8610-c3ff318b117b
# ╟─3e19d251-91f6-4383-bfc8-ffe816570f42
# ╟─ce4f2fdd-16b1-46e8-88a4-a952896b6df8
# ╟─dd1561e2-233f-425a-832f-130b49f0bf0b
# ╟─05ec0ea2-cfec-4e91-a46a-68bbdaefd562
# ╠═e480137d-b6d9-4e18-92f0-640292bbb5f0
