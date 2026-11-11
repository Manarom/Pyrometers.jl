### A Pluto.jl notebook ###
# v1.1.0

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
	notebook_dir = @__DIR__()
	src_dir = joinpath(abspath(joinpath(notebook_dir,"..")),"src")
	import Pkg 
	Pkg.activate(notebook_dir)
	#Pkg.resolve()
	using Revise
	using Pyrometers
	using PlutoUI
	using DelimitedFiles
	using Interpolations
	using Plots , PrettyTables
	using JLD2
	using NumericalIntegration
	using OrderedCollections
end;

# ╔═╡ 8911823a-2210-4328-b9a5-aa88af2146c5
using StatsBase

# ╔═╡ b1e73e8f-8eb4-4688-86a8-b787d6b15454
using QuadGK

# ╔═╡ 4f1397e1-018d-4965-b3e1-78ccac62b54c
using BenchmarkTools

# ╔═╡ 30743a02-c643-4bdc-837e-b97299f9520a
md"""
#  `Pyrometry.jl`   

######  **`Pyrometry.jl`** содержит набор методов и типов для воздания виртуальных пирометров. Виртуальный пирометр может "измерять" температуру по загружаемой интенсивности излучения, а также подстраивать излучательную способность пирометра под заданную "измеренную" температуру.
____________________
## Установка пакета

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
3) Copy the entire GitHub [repository](https://github.com/Manarom/RadiationPyrometry.jl.git) to your local folder
4) Open this notebook file located at  `project_folder\notebooks\radiation_pyrometry.jl` in `Pluto` by providing the full path to the *"Open a notebook"* text field on `Pluto`'s starting page.

"""

# ╔═╡ ba23c985-74c4-41f3-8bc4-f7287e30e47f
#using Revise,StaticArrays,OrderedCollections,Optimization,OptimizationOptimJL,LaTeXStrings,NumericalIntegration,Interpolations,Plots,PlutoUI,DelimitedFiles , ForwardDiff , PlutoPlotly , PrettyTables, JLD2 , Roots , QuadGK , ADTypes , ForwardDiff , StatsBase

# ╔═╡ 89a11dcd-b3b5-4349-930d-a66ad74e8fa2
const Planck = Pyrometers.Planck

# ╔═╡ 171409eb-22b5-4bc5-a8e2-eac0932a24f3
PlutoUI.TableOfContents(indent=true, depth=4, aside=true)

# ╔═╡ 04f9f75e-34cc-4921-971e-6958ed960518
begin 
	plot_common_args = (grid = true, gridlinewidth=3, gridstyle = :dot,minorgrid=true, box = :on, linewidth = 3 , fmt = :png)
	data_folder = joinpath(notebook_dir , "data")
	data_names = readdir(data_folder)
end

# ╔═╡ 7d37ae06-c81b-4aa1-a05c-e03b83e395d7
data_folder

# ╔═╡ 7cfee16a-4c47-4ded-9b5a-4a21e59dd9ad
begin
html"""
    <style>
        /* All input[type="range"] sliders */
        input[type="range"] {
            width: 100%;          /* or e.g. 300px */
            min-width: 600px;     /* optional cap */
        }
        input[type="text"] {
            width: 100%;          /* or e.g. 300px */
            max-width: 600px;     /* optional cap */
        }
        /* Or target specific sliders by class */
        .my-slider {
            width: 400px;
        }
    </style>
"""
end

# ╔═╡ d5ee3913-66be-47d7-a755-699ba64b4f98
md"""
## Введение

Этот ноутбук содержит три примера использования пакета **Pyrometers.jl**. Цель этого небольшого пакета — создание виртуального пирометра, который можно использовать для расчета коэффициента излучения реального прибора. Это позволяет корректировать измеренную температуру так, чтобы она соответствовала фактической температуре нагретого объекта.

"""

# ╔═╡ d442014a-20e6-4be4-ac7f-f13de329dec5
md"""
## I. Тепловое излучение: `абсолютно черного тела` и `реальной поверхности` 
_______________________

Все нагретые тела испускают тепловое излучение. Согласно закону Планка, спектр идеального излучателя (называемого *абсолютно черным телом*) определяется исключительно его температурой. Спектральная интенсивность излучения абсолютно черного тела рассчитывается следующим образом: \


``I_{blackbody}(\lambda , T) =  \frac{C_1}{\lambda ^5} \cdot \frac{1}{e^{\frac{C_2}{\lambda T } } - 1}``, \


где ``C_1`` = $(Planck.C₁), ``Вт \cdot мкм/м² \cdot ср`` и ``C_2`` = $(Planck.C₂), ``мкм \cdot К``, ``\lambda`` — длина волны в ``мкм``, ``T`` — температура в Кельвинах.

Интенсивность теплового излучения реальной поверхности ниже, чем у абсолютно черного тела. Доля интенсивности излучения, испускаемого реальной поверхностью по сравнению с абсолютно черным телом, характеризуется спектральной направленной излучательной способностью ``\epsilon (\lambda, T,\vec{\Omega})``:

``I_{real\ \ surface}(\lambda , T, \vec{\Omega} ) = \epsilon (\lambda, T,\vec{\Omega}) \cdot  I_{blackbody}(\lambda , T)``, \


здесь ``\vec{\Omega}`` обозначает направление. \


Интересно, что, в отличие от абсолютно черного тела, тепловое излучение реальной поверхности (в общем случае) зависит от направления излучения. Таким образом, наиболее общей характеристикой теплового излучения реальной поверхности является именно `спектральная направленная излучательная способность`.

В пакете [PlanckFunctions.jl](https://manarom.github.io/PlanckFunctions.jl) содержится ряд функций для расчета спектров теплового излучения абсолютно черного тела (а также различных производных, интегралов и т.д.).
На следующем графике показано влияние спектральной степени черноты на интенсивность теплового излучения реальной поверхности.
"""

# ╔═╡ 27b3c586-9eb0-4a51-b9ca-a9c0379fccdf
@bind  λ_BB PlutoUI.combine() do Child
	md"""
	Blackbody spectral range, ``\mu m`` : \
	``\lambda_{left}`` = $(
		Child(Slider(0.1:0.1:50,default=0.1,show_value = true))
	)   -- 
	 $(
		Child(Slider(0.1:0.1:150,default=15.0,show_value = true))
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

# ╔═╡ 51de94a8-60ae-4099-97f8-f20582c4bc48
md""" 
Что построить : $(@bind what_to_plot Select([:i,:di,:d2i]))
"""

# ╔═╡ 72095f49-f6c0-41a8-be2f-323a0bacfc51
md""" Материал поверхности: $(@bind material_selection Select(data_names , default = "zrb2.txt"))"""

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

# ╔═╡ e64ba0f6-48a6-4e87-a347-e824f350ef4a
begin 
    # 1. Задаем заголовки колонок
    atmospheric_windows_headers = [
        "Окно прозрачности", 
        "Название диапазона", 
        "Спектр (мкм)", 
        "Газы-поглотители", 
        "Тип физического приемника", 
        "Диапазон температур",
        "Примеры моделей пирометров"
    ]
    # 2. Заполняем матрицу данных
    atmospheric_windows_data  =[
        "Ближнее ИК (SWIR)"   "Коротковолновое ИК"  "0.95 - 1.35\n1.50 - 1.80" "H2O, CO2"              "Кремниевый фотодиод (Si)\nИндий-галлий-арсенид (InGaAs)" "от +250 °C до +3000 °C+" "• Modline 51, 61, 71\n• Optris CT / CTlaser 1M\n• Lumasense Impac IS 5, IS 6, IS 140";
        "Ближнее ИК (SWIR-2)" "Средне-короткое ИК"  "2.00 - 2.40"              "H2O, CO2"              "Расширенный InGaAs\nСульфид свинца (PbS)"               "от +50 °C до +1800 °C"   "• Modline 52, 5G\n• Optris CT / CTlaser 2M, 3M\n• Lumasense Impac IGA 5, IGA 6";
        "Среднее ИК (MWIR)"   "Средневолновое ИК"   "3.00 - 5.00"              "H2O, CO2\n(провал 4.2-4.4)" "Селенид свинца (PbSe)\nАнтимонид индия (InSb)"       "от +50 °C до +1200 °C"   "• Modline 43 (пластик)\n• Optris CT G5, CTlaser P7\n• Lumasense Impac IPE 140, IN 5/5";
        "Длинноволновое (LWIR)" "Главное тепловое"  "8.00 - 14.00"             "H2O, O3\n(провал 9.6)" "Термобатареи (Thermopiles)\nМикроболометры"              "от -50 °C до +900 °C"    "• Modline 44, Raytek GP\n• Optris CS LT, CSmicro LT, CT LT\n• Lumasense Impac IN 5, IN 300"
    ]

    # 3. Отрисовка таблицы в HTML для красивого отображения внутри Pluto
    pretty_table(
        HTML, 
        atmospheric_windows_data; 
        column_labels = atmospheric_windows_headers,
        alignment = :l, # Выравнивание по левому краю для читаемости текста
        line_breaks = true, # Разрешаем перенос строк (\n) внутри ячеек
        #tf = PrettyTables.tf_html_minimalist # Аккуратный минималистичный стиль без лишней пестроты
    )
end

# ╔═╡ 2ad3ec82-54a2-49ac-94ef-579f808dfb1a
begin 
	rt_emissivity_data = readdlm(joinpath(data_folder , material_selection)) # loading file, actually this file is two-column data with no headers, but this is ok
	rt_emissivity_interpolation = linear_interpolation(rt_emissivity_data[:,1],rt_emissivity_data[:,2],extrapolation_bc = Interpolations.Flat())
end;

# ╔═╡ 4d6337aa-cfc7-4154-a395-5aa53e23d01a
begin 

	
	T₁ = T_BB[1] + Planck.Tₖ # converting to Kelvins
	T₂ = T_BB[2] + Planck.Tₖ # converting to Kelvins
	λ_bb = collect(range(λ_BB...,length=500))
	ll = 1e-3
	if what_to_plot == :i 
		f_to_plot = Planck.ibb
	elseif what_to_plot == :di
		f_to_plot = Planck.:∇ₜibb
	else 
		f_to_plot = Planck.:∇²ₜibb
	end
		
	
	ibb1 = [  ifelse(i>ll, i , ll) for i in  f_to_plot.(λ_bb,T₁) ]
	ibb2 =  [  ifelse(i>ll, i , ll) for i in f_to_plot.(λ_bb,T₂) ]
	ireal1 = ibb1.*rt_emissivity_interpolation(λ_bb)
	ireal2 = ibb2.*rt_emissivity_interpolation(λ_bb)
	mmx = max([maximum(i) for i in (ibb1,ibb2)]... )
	fls = min(minimum(ireal1),minimum(ireal2))
	
	p_bb = Plots.plot(λ_bb, ibb1 ,label="АЧТ T=$(T₁),K", xscale=scales_BB[1],yscale =scales_BB[2], linewidth = 3)#,fillrange=fls, fillalpha=0.3; plot_common_args... )
	
	plot!(λ_bb , ibb2 , label="АЧТ T=$(T₂),K", xscale=scales_BB[1],yscale =scales_BB[2] , linewidth=3) #, fillrange=fls, fillalpha=0.3 )
	
	Plots.plot!(λ_bb,ireal1,label="$(replace(material_selection , ".txt"=>"") ) : T=$(T₁),K", xscale=scales_BB[1],yscale =scales_BB[2],fillrange=fls, fillalpha=0.3; plot_common_args...)	
	
	Plots.plot!(λ_bb,ireal2,label="$(replace(material_selection , ".txt"=>"") ) : T=$(T₂),K", xscale=scales_BB[1],yscale =scales_BB[2],fillrange=fls, fillalpha=0.2; plot_common_args...  , legend = :best , fmt=:png)
	xlabel!("Длина волны, μm")
	ylabel!("Спектральная интенсивность, W/m²⋅sr⋅μm")
	ylims!((1e-3 , mmx*1.1 ))
end

# ╔═╡ f5a2937c-ce1c-4c03-8e86-8db363bd21c9
begin 
	rel_error_I = 1e-2
	p_error = Plots.plot()
	for T₁ in 273.15 : 200 : 2273.15
	#rel_error_T = rel_error_I ./ (T₁ .* Planck.∇ₜibb.(λ_bb, T₁) ./ Planck.ibb.(λ_bb, T₁))
	rel_error_T =[ ifelse(i>ll, i , ll) for i in  Planck.∇ₜibb.(λ_bb, T₁)]
	Plots.plot!(p_error , λ_bb , rel_error_T;plot_common_args... , label = "T=$(round(T₁ , digits=1)) K" , linestyle = :auto , yscale = :log10)
	end
	p_error
	xlabel!(p_error,"Длина волны, мкм")
	ylabel!(p_error,"db(λ,T)/dT, $(Planck.units(Planck.∇ₜibb))")
end

# ╔═╡ 1e429893-8aaa-4e95-aff8-4eab826a5c4f
fls

# ╔═╡ 03d76e64-ebf4-432b-b9be-d4cb26275f55
begin 
	e_real_BB = rt_emissivity_interpolation(λ_bb) # this spectral emissivity was measured up to 18 μm, thus for higher wavelengths it uses flat extrapolation
	Plots.plot(λ_bb, e_real_BB, label=nothing, linewidth=3.0; plot_common_args...)
	title!("Surface spectral emissivity of $(material_selection)")
	xlabel!("Wavelength, μm");ylabel!("Spectral emissivity (ϵ)")
end

# ╔═╡ 8a066ee5-80e9-462f-9a61-15851468aa63
md"""
## II. Пирометрия частичного излучения
_______________________

Поскольку энергия теплового излучения `абсолютно черного тела` сильно зависит от температуры, эту величину можно использовать для измерения температуры реальной поверхности. В этом заключается основная идея пирометрии частичного излучения: **измерить интенсивность, чтобы получить температуру**. Так как интенсивность является величиной направленной, пирометру требуется коллимирующая оптика (телескоп). Степень черноты реальных поверхностей часто существенно меняется в зависимости от длины волны; в то же время пирометры частичного излучения предполагают постоянную степень черноты (так называемое «серое» приближение в полосе пропускания). Таким образом, для промышленных целей полезно иметь несколько пирометров, каждый из которых работает в относительно узком спектральном диапазоне. В спектральном диапазоне пирометра частичного излучения степень черноты не должна значительно изменяться, чтобы допущение о её постоянстве оставалось актуальным.

Пакет [Pyrometers.jl](https://manarom.github.io/Pyrometers.jl) предоставляет ряд функций для работы с `виртуальными` пирометрами частичного излучения.

"""

# ╔═╡ b7fac177-c211-4635-992f-e6473be7bdae
md"""
Словарь **Pyrometers.DefaultPyrometersTypes** содержит названия стандартных типов пирометров вместе с их спектральными диапазонами. Пользовательский пирометр можно создать, указав его тип (имя), длину волны или спектральный диапазон, а также рабочую степень черноты:

```julia
# создание пользовательских объектов-пирометров
p_custom = Pyrometers.Pyrometer((2.5, 3.7) , type = :Custom, ϵ = 0.65)
```
"""

# ╔═╡ 69f9052b-d925-49fb-8eb9-089671e64339
md" ### II.I Типовые пирометры частичного излучения "

# ╔═╡ e2a9aa39-2490-4681-89d3-a01f058f6feb
pretty_table(HTML, [String(p[1]) => p[2] for p in Pyrometers.DefaultPyrometersTypes],top_left_string ="Таблица типовых пирометров в пакете Pyrometers.jl и соотвествующие диапазоны длин волн", wrap_table_in_div=true )


# ╔═╡ 02eee968-ff43-4b82-8d68-efede1a220dd
md"""
После создания объекта **Pyrometer** его можно использовать для «измерения» температуры, а именно — для преобразования интенсивности излучения реальной поверхности в значение температуры в соответствии со спектральным диапазоном прибора. Это делается с помощью функции ``Pyrometers.measure(pyr, measured_intensity)`` (интенсивность должна быть указана в соответствующих единицах [Вт/м²⋅ср⋅мкм]) или просто `pyr(measured_intensity)`, так как объект - пирометр является вызываемым. На практике часто бывает известна истинная температура поверхности в какой-то момент; в этом случае степень черноты виртуального пирометра можно скорректировать так, чтобы «измеренная» температура совпала с реальной. Это можно сделать с помощью функции:

```julia
Pyrometers.fit_ϵ!(p::Pyrometer, Tmeasured::Float64, Treal::Float64) 
```

"""

# ╔═╡ 6c2d1616-4992-4123-a356-0545cf81ce08
md" ### II.II Измерение пирометром без задания излучательной способности "

# ╔═╡ 9b08b767-7e8f-4483-9f2f-226022ce10e4
md"""
Если излучательная способность пирометра частичного излучения задана неверно, измеренная температура также будет неточной. Интересно посмотреть, как различные пирометры (с установленной единичной степенью черноты) «измеряют» температуру реальной поверхности. На следующем графике показана интенсивность теплового излучения реальной поверхности и рабочие диапазоны нескольких распространенных типов пирометров. В легенде указана их «измеренная» температура. Значение «измеренной» температуры для каждого типа пирометра получено путем сопоставления мощности излучения абсолютно черного тела с мощностью излучения реальной поверхности, проинтегрированных в пределах рабочего спектрального диапазона прибора.
"""

# ╔═╡ 712828a7-fb54-42e6-95fc-233243190f59
md"Real surface temperature $(@bind T_pyr Slider(range(10,3000,1000),default=1500,show_value=true) ) "

# ╔═╡ c69acbf6-94fb-4ac3-8d56-d1f9dda11440
begin
	λ_pyr = collect(range(0.1,18,1000))
	pyrometers_vector = sort(Pyrometers.produce_pyrometers())# returns a vector of all default pyrometers
	N = length(pyrometers_vector) + 1
	# calculationg the real surface thremal radiation spectrum
	real_i = Planck.ibb.(λ_pyr,T_pyr).*rt_emissivity_interpolation(λ_pyr)
	
	data_legend = Matrix{String}(undef,N,1)
	data_legend[1] = "T_{real} =$(round(T_pyr))"
	# poltting surface thremal emission spectrum
	plot_pyrometers = Plots.plot(λ_pyr,real_i,xscale=scales_BB[1],yscale =scales_BB[2],fillrange=0, fillalpha=0.3,dpi=600,label = data_legend[1],legend_background_color=:white,legend_foreground_color = :black,legend_position=:right; plot_common_args...)

	real_i_interp = linear_interpolation(λ_pyr,real_i)
	xlabel!("Длина волны , μm")
	ylabel!("Интенсивность излучения")
	#ylabel!()
	max_val = maximum(real_i)
	t_em_unity = Vector{Float64}(undef,length(pyrometers_vector))
	t_em_acttual = Vector{Float64}(undef,length(pyrometers_vector))
	for j in 1:length(pyrometers_vector)# iterating over virtual pyrometers vector
		# the following function checks if current pyrometer is narrow band
		ppp = pyrometers_vector[j]
		is_two_wavelength_pyrometer = Pyrometers.is_spectral_band(ppp)
		# 
		l_cur = is_two_wavelength_pyrometer ? ppp.λ : [ppp.λ[]-0.2,ppp.λ[]+0.2 ]
		λ_pyr_interp = collect(range(l_cur...,length=30))
		# calculating the measured by the pyrometer value 
		measure_intensity =  is_two_wavelength_pyrometer ? NumericalIntegration.integrate(λ_pyr_interp,real_i_interp(λ_pyr_interp)) : real_i_interp(ppp.λ[1])
		# setting emissivity to one
		Pyrometers.set_emissivity!(ppp,1.0)
		measured_temp = Pyrometers.measure(ppp,measure_intensity,T_starting = T_pyr)
		# temperature measured by the current pyrometer 
		data_legend[j+1] = "$(ppp.type): T="*string(round(measured_temp))
		# plotting current pyrometer spectral range
		region_flag = 
		plot!(l_cur,[max_val,max_val],fillrange=0, fillalpha=0.5,label=data_legend[j+1]; plot_common_args...)
		# remember the value of temperature with unit emissivity
		t_em_unity[j] = measured_temp 
		# calculating the averaged gray-band emissivity
		Pyrometers.fit_ϵ!(ppp,measured_temp,T_pyr)
		# calcaulting the averaged emissivity wthin the pyrometers spectral band (or at fixed wavelength)
		 t_em_acttual[j] =  round(Pyrometers.measure(ppp,measure_intensity,T_starting= T_pyr))
		 
	end
	plot!(twinx(),λ_pyr, rt_emissivity_interpolation(λ_pyr),linewidth=4,linecolor=:red,label=nothing,alpha=0.3,ylabel ="Излучательная способность поверхности"; plot_common_args... )
	
	
	# emissivities table 
	data = Matrix{Any}(undef,length(pyrometers_vector),5)
	e_grey = [p.ϵ[] for p in pyrometers_vector]
 	data[:,3:end] .= hcat(t_em_unity,e_grey, t_em_acttual)
	data[:,1] .= [p.type for p in pyrometers_vector]
	data[:,2] .= [string(p.λ) for p in pyrometers_vector]

end;

# ╔═╡ 7071a6f4-e296-4e53-8e6f-24f1f038c1a5
plot_pyrometers

# ╔═╡ f763d449-2a7a-4008-a183-823a774bc25e
#savefig(plot_pyrometers,joinpath(notebook_dir,"Pyrometers.png"));

# ╔═╡ 667f7c30-56e0-461f-b35b-c924007eb9f2
md"""
Для данного конкретного материала показания пирометра типа **F** наиболее близки к реальной температуре, так как степень черноты реальной поверхности в спектральном диапазоне этого прибора ближе к единице. Все результаты сведены в следующую таблицу.
"""

# ╔═╡ a861d56f-f6c9-4754-b7a9-ed63713f1f2f
	pyr_table = pretty_table(HTML,data, column_labels= ["тип","λ ,μm","T₀ (ϵ=1),K","grey-ϵ", "T₁ (grey-ϵ),K"],top_left_string ="Таблица температур, `измеренных` различными пирометрами при излучательной способности заданной равной единице (T₀), рассчитанная излучательная способность, которая должна быть задана на пирометре (grey-ϵ) и показания пирометра после задания этой излучательной способности (T₁), истинная температура поверхности Tᵣ=$(T_pyr)" , wrap_table_in_div=true)

# ╔═╡ d08ec8f2-7689-4043-9f28-da06ab0124b9
md"""
## III. Излучательная способность эталонного источника абсолютно черного тела
_______________________

Реальный пирометр не измеряет интенсивность излучения напрямую из-за спектральной зависимости его чувствительности и множества других факторов; поэтому он требует калибровки. Чтобы откалибровать пирометр, его помещают перед эталонным источником, тепловое излучение которого максимально приближено к излучению абсолютно черного тела. Обычно таким эталоном служит специально сконструированная печь с нагреваемой полостью и точным контроллером температуры.

В процессе калибровки оператор пирометра измеряет температуру эталонного источника. В идеале измеренная температура должна совпадать с температурой, установленной на контроллере эталона, однако на практике эти значения редко совпадают точно. Это расхождение возникает из-за неидеальной природы эталонного источника. Спектральная излучательная способность реального эталона не равна единице и варьируется в зависимости как от длины волны, так и от температуры. Чтобы учесть это отклонение, источник поставляется вместе с калибровочной таблицей, которая выглядит следующим образом:


"""

# ╔═╡ 3ddc9d4e-ff96-4173-899b-3eae0facb830
md" ### - Калибровочная таблица эталонного излучателя "

# ╔═╡ c5ac80ee-8143-4c28-bbff-2ac761c71fac
begin
	bb_calibration_table_data = readdlm("BBethalon")
	table_header =["Tref";]
	all_types = [getfield(p,:type) for p in pyrometers_vector]
	table_header = vcat(table_header,all_types)
	pr_tbl = pretty_table(HTML,bb_calibration_table_data,column_labels= table_header,top_left_string = "Пример калибровочной таблицы для эталона АЧТ (все температуры  - в градусах Цельсия)",  wrap_table_in_div=true)
	bb_calibration_table_data .+= Planck.Tₖ # converting table data to 
	ref_T = @view bb_calibration_table_data[:,1] # reference source temperature
	pr_tbl
end

# ╔═╡ 0c9fe7b1-374c-4fb8-9cfe-9337389713bf
md"""
### - Определение спектра излучательной способности по таблице калибровочных температур

Первая строка в таблице представляет температуры эталона, в то время как последующие столбцы показывают температуры, измеренные различными пирометрами. Каждый столбец помечен в соответствии с типом пирометра. Следует отметить, что значения температур в таблице меняются как в зависимости от типа пирометра (в каждой строке), так и от самих значений температуры (в каждом столбце). Это указывает на то, что спектральная излучательная способность эталонного источника зависит как от длины волны (соответствующей типу пирометра), так и от температуры. Поскольку различные типы пирометров в совокупности охватывают широкий спектральный диапазон от 2 до 14 мкм, данные измеренных температур могут быть использованы для определения спектральной излучательной способности эталона и её температурной зависимости на основе приведенной выше таблицы калибровки. Для этого излучательная способность пирометра подбирается таким образом, чтобы показываемая им температура в точности соотвествовала температуре АЧТ.

На следующем рисунке показана спектральная излучательная способность эталона абсолютно черного тела, рассчитанная по данным таблицы калибровки температуры.

"""


# ╔═╡ bc2d93ae-6c30-462c-96e0-30fdb84d7c63
md"""
Три температуры АЧТ: 

$(@bind T_ref1  Slider(ref_T,show_value=true,default = ref_T[1])) 

$(@bind T_ref2  Slider(ref_T,show_value=true,default = ref_T[1])) 

$(@bind T_ref3  Slider(ref_T,show_value=true,default = ref_T[end])) 

"""

# ╔═╡ d3199b6e-9779-4def-b701-fe85d1035045
begin 
	#wavlength = RadiationPyrometers.wlength
	full_wavelengths_range = Pyrometers.full_wavelength_range(pyrometers_vector)
	jj = indexin(T_ref1,ref_T)[]
	foreach(pyrometers_vector) do p
			Pyrometers.set_emissivity!(p , 1.0)
	end
	ej = Pyrometers.fit_ϵ_wavelength!(pyrometers_vector,T_ref1,bb_calibration_table_data[jj,2:end])
	pppp = Plots.plot(full_wavelengths_range, ej, title="Spectral emissivity of the blackbody reference",label ="T = $(ref_T[jj]) , K",grid=true; plot_common_args...)
	for T_reference in (T_ref2 , T_ref3)
		foreach(pyrometers_vector) do p
			Pyrometers.set_emissivity!(p , 1.0)
		end
		global jj = indexin(T_reference,ref_T)[]
		global ej = Pyrometers.fit_ϵ_wavelength!(pyrometers_vector , ref_T[jj],bb_calibration_table_data[jj,2:end])
		Plots.plot!(pppp,full_wavelengths_range , ej, label=" T = $(ref_T[jj]) , K"; plot_common_args...)
	end
	xlabel!(pppp,"Wavelength, μm")
	ylabel!(pppp,"Emissivity")
	pppp
end

# ╔═╡ 6aa82819-6f5a-43af-ac03-4775fce51564
md"""
Из рисунка выше видно, что излучательная способность эталона изменяется с температурой слабо и наиболее близка к единице в области блин волн меньше `3` и больше `7` μm  
"""

# ╔═╡ 72947d97-0a97-4064-a2fe-08d19dec0f0e
md"""
 ## IV Измерение температуры поверхности при наличии "паразитной" засветки
"""

# ╔═╡ 854f8c06-aac6-464f-899c-41ca706b259a
md"""

На практике возникают ситуации, когда необходимо производить измерения температуры поверхности в условиях внешней застветки. Например, это может быть актуально при нагреве поверхности высокоинтенсивным источником излучения. 

"""

# ╔═╡ 9ce196c1-8915-46da-9aba-f13d7655959a
md"""
 Теперь излучение пропорционально спектральной интенсивности двух источников абсолютно черного тела (АЧТ) с разными температурами. Излучение паразитного источника отражается от поверхности образца с коэффициентом отражения **``\rho = 1 - \epsilon``** и попадает в пирометр.  Предполагается, что паразитное излучение напрямую в пирометр не попадает. 

Таким образом, интенсивность регистрируемая пирометром может быть выражена следующей формулой:

``I_{meas} = \epsilon(λ)\mathcal{b}(\lambda , T_1) + (1 - \epsilon(λ))\mathcal{b}(\lambda , T_2)``

Здесь, ``\epsilon(λ)`` - это излучательная способность поверхности 

"""

# ╔═╡ 86af6afc-b28a-4e84-952a-bd29710374f8
λ2 = collect(range(0.8,17.0,1000));

# ╔═╡ 620dce38-97ce-495f-9b23-1b8290cbd973
begin 
	to_names = Dict(data_names[1] => "Cr2O3" , data_names[2] => "Black"  , data_names[3] => "DKS#8", data_names[4] => "кварцевая керамика" ,data_names[5] => "РСНК"   , data_names[6] => "ZrB2" , "fixed" => "fixed")
end;

# ╔═╡ b15bd3b0-0971-4b22-a109-0ed56afda615
begin 
	all_emissivities = Plots.plot(;plot_common_args...)
	all_emissivities_data = OrderedDict{String,AbstractExtrapolation}()
	λ_plot = λ2[1:20:end]
	for et in  data_names
		_data = readdlm(joinpath(data_folder , et))
		eint = linear_interpolation(_data[:,1] , _data[:,2] , extrapolation_bc=Line())
		ϵ = eint.(λ_plot)
		Plots.plot!(all_emissivities,λ_plot , ϵ ; plot_common_args... , label = to_names[et] , linestyle = :auto , markershape = :auto , markerstrokecolor = :auto , markersize = 8 )
		push!(all_emissivities_data , et => eint)
	end
	all_emissivities
	xlabel!(all_emissivities , "Длина волны, мкм")
	ylabel!(all_emissivities , "Излучательная способность")
end

# ╔═╡ 31a59906-e0d8-4ae8-b1ca-34b9f01062c8
all_emissivities_data

# ╔═╡ 0cb1cfa0-6f2d-49bc-ac35-bb76688bb510


# ╔═╡ ac15bdd0-7679-4511-b0be-4885f2d366bd


# ╔═╡ d67ef942-be2d-4fad-91ba-eff9eda9d8e3
begin 
	p_in = Pyrometers.SpectralBandPyrometer(8.0 , 9.7 , type = :IN59)
	p_thermo = Pyrometers.SpectralBandPyrometer(8.0 , 14.0 , type = :UV640)
	special_pyrometers = (p_in , p_thermo)
end

# ╔═╡ da0a4a00-692a-4d29-97cb-74366e8374c9
begin 
	e_int_plot = (Plots.plot(;plot_common_args...) , Plots.plot(;plot_common_args...) )
	T_integra = range(273,1973,20)
	
	for (pyr_number , pyr) in enumerate(special_pyrometers)
		for (k , d) in all_emissivities_data
			_em_wrapper = Pyrometers.IsothermalSpectralQuantity(d)
			_e_int = fill(0.0 , length(T_integra))
			for (i , t) in enumerate(T_integra)
				_e_int[i] = first(Pyrometers.integral_emissivity(pyr , _em_wrapper , t))
			end
			plot!(e_int_plot[pyr_number] , T_integra .- 273.15, _e_int;
				  plot_common_args... , xlabel = "Температура, ᵒC" , ylabel = "Интегральная ИС" , label ="$(to_names[k]) : <ϵ> = $(round(StatsBase.mean(_e_int) , digits = 2))" , legend = :bottomright ,  linestyle = :auto , markershape = :auto , markerstrokecolor = :auto , markersize = 8 )
			title!(e_int_plot[pyr_number] , "$(Pyrometers.shorthand(pyr))")
		end
	end
	#title!()
end

# ╔═╡ 90e15591-d885-48f3-bc46-cb9385c9ed13
ylims!(e_int_plot[1] , (0.8 , 1.01))

# ╔═╡ afb61470-2455-4742-b81f-63baf5a4c3ce
ylims!(e_int_plot[2] , (0.65 , 1.01))

# ╔═╡ 8cef05a1-2974-4c38-b73e-fa706f347fcc
md" Show wavelength range: $(@bind λ_show RangeSlider(range(extrema(λ2)... , 1000)))"

# ╔═╡ 36ba2396-bb5e-4d22-a58e-9ab27cd18b2d
md" Излучательная способность поверхности , ϵ = $(@bind ϵ_fixed  Slider(1e-4:1e-4:1.0 , show_value = true , default = 0.5))"

# ╔═╡ 91bbd553-4e4a-431d-9d54-b0f4882fd426
md"""
Температура образца  ``T_1 `` = $(@bind T1 Slider(300.0:1e-2:3000 , show_value = true , default = 1273.15)), K

"""

# ╔═╡ 1c2bd109-8023-49a5-87ba-32d755b201dd
@bind selected_standard_pyros MultiSelect(collect(keys(Pyrometers.DefaultPyrometersTypes)) , default=collect(keys(Pyrometers.DefaultPyrometersTypes)))

# ╔═╡ ad622962-9cd8-432c-a988-ec8b3d676077
md" Использовать излучение объекта в виде интерполятора $(@bind use_external_lamp CheckBox(true))"

# ╔═╡ 15f1519b-d924-4fa9-b212-eebba75c544a
if !use_external_lamp 
md"""
Температура паразитного источника ``T_2`` = $(@bind T2 Slider(300.0:1e-2:5000 , show_value = true , default = 3280.0)), K
"""
else
md"""
Напряжение на лампе ``U`` = $(@bind T2 Slider(100.0:1:240 , show_value = true , default = 220.0)), В
"""
end

# ╔═╡ c918b685-1b6e-4ce3-a1e6-501c6b4ffe0e
md" Interpolator location $(@bind interpolator_folder TextField(120,default = @__DIR__()))"

# ╔═╡ 96feae26-99e4-45c9-869b-410ec7965b63
interpolators_in_folder = collect(f for f in readdir(interpolator_folder) if contains(f , ".jld2")  );

# ╔═╡ 1811e43c-f7db-47b1-9b83-bb38455d7db3
#pyrometers_vector2 = deepcopy(pyrometers_vector);

# ╔═╡ 8ece7476-4257-49f8-a400-50bcc347908d
use_external_lamp && !isempty(interpolators_in_folder) && md""" Select interpolator file $(
@bind interpolator_file_name Select(
interpolators_in_folder
)
)"""

# ╔═╡ c5a00df1-ff02-4564-9838-e678f8ea7389
if use_external_lamp && !isempty(interpolators_in_folder)
	lamp_interpolator =JLD2.load(joinpath(interpolator_folder , interpolator_file_name))["lamp_interpolator"]
	voltage_range = [100.0 , 140.0 , 180.0 , 220.0 , 240.0]
end

# ╔═╡ 8b441a27-004a-4468-b89f-f3a6cc6f9bbf
begin #integral absorptivity with respect to the external radiation 
	v = 100:20:240
	A = OrderedDict{String , Vector{Float64}}()
	for (k , e_int) in all_emissivities_data
		_A_int = Vector{Float64}(undef , length(v))
		_e_surf = Pyrometers.IsothermalSpectralQuantity(e_int)
		for (i , vi) in enumerate(v) 
				_lamp = Pyrometers.IsothermalSpectralQuantity(l->lamp_interpolator(l,vi))
				denom = quadgk(_lamp , 0.4 , 15.0)
				enum_iso = Pyrometers.fix_temperature(_e_surf * _lamp , 1000.0)
				_enum  = quadgk(enum_iso , 0.4 , 15.0)
				_A_int[i] = first(_enum)/ first(denom)
		end
		A[k] = _A_int
	end


end

# ╔═╡ 6edddd01-915c-45bb-b3ef-155f4dcd9d8c
first(A)

# ╔═╡ bfce276c-389c-405c-a552-979ed38133a2
begin 
	p_A_int = Plots.plot(;plot_common_args...)
	for (k,d) in A 
		_a_avg = round(mean(d) , digits = 2)
		nm = to_names[k]
		Plots.plot!(p_A_int , v , d; plot_common_args... , label = nm*"($_a_avg)" , markersize=8 , markershape = :auto  , legend_position = :top)
	end
	p_A_int
	ylims!(p_A_int , 0 , 1.3)
	xlabel!(p_A_int , "Напряжение  на лампе, В")
	ylabel!(p_A_int , "Поглощательная способность")
end

# ╔═╡ 0f0baa19-aa64-4daa-a501-861e22cdafde
lamp_interpolator(2.1, 120)

# ╔═╡ f3373768-d598-4d1e-b01b-1929e233f36a
lamp_interpolator(2.1, 220)

# ╔═╡ b27eff45-93e2-4aeb-a7ef-88132dd7e9e8
ll2 = Pyrometers.IsothermalSpectralQuantity(l->lamp_interpolator(l,230))

# ╔═╡ 4cc48476-fb0c-4a30-a93d-2e6caa53e232
first(Pyrometers.quadgk(Pyrometers.fix_temperature(0.7*ll2 , 1400) , 0.1 , 15.0))/first(Pyrometers.quadgk(ll2 , 0.1 , 15.0))

# ╔═╡ 54339700-71fd-48bf-a2ef-0c3267b9d81b
md" **Пересчитать матрицу ошибки измерений $(@bind is_recalculate CheckBox(false))**"

# ╔═╡ 7f76bc22-77c3-4eb3-9d49-588e653df2e7
md"""

###### Расчет ошибки определения температуры поверхности для выбранного типа пирометра для различных значений температур паразитного излучения и поверхности   



"""

# ╔═╡ a0197a9a-34bf-4a3e-af8a-c23ea777f482
md" **Пользовательский диапазон длин волн** : $(@bind use_custom CheckBox(default = false))"

# ╔═╡ a90bde8a-222b-4b93-892b-9fa4cfc43372
md" **Пользовательский пирометр** : $(@bind use_custom_pyros CheckBox(default = true))"

# ╔═╡ d86d7e8d-94e7-4d1b-b4d8-df52e338935d
#md" ИСпользовать заданную ИС $()

# ╔═╡ 7793976e-6714-4bc1-9d97-412dd2a67480
@bind custom_waves PlutoUI.combine() do Child

md"""
	
Диапазон длин волн пользовательского пирометра:

λleft = $(Child("left", NumberField(0.1:1e-3:20 , default = 8.0)))

λright = $(Child("right", NumberField(0.1:1e-3:20 , default = 9.7)))

"""
end

# ╔═╡ c01e05c7-2b29-4f57-96ee-9a4eebef9025
md"""

User defined pyrometer $(@bind custom_selection Select( [p_in=>"IN59" , p_thermo=>"UV640"]  ))
"""

# ╔═╡ f181980f-bf72-4468-8daa-9461c6c901e0
md"""
Можно выбрать различные излучательные способноти : $(@bind  emissivity_type Select(vcat(["fixed"] , data_names), default = "fixed"))

"""

# ╔═╡ efc35420-d0e6-4795-94b6-d43289b4de44
begin 
	if emissivity_type == "fixed"
		ϵ = ϵ_fixed
		eint = Returns(ϵ_fixed)
	else
		_data = readdlm(joinpath(data_folder , emissivity_type))
		eint = linear_interpolation(_data[:,1] , _data[:,2] , extrapolation_bc=Line())
		ϵ = eint.(λ2)
		 
	end
end;

# ╔═╡ a7ff8a2d-a81d-4474-b27f-565de2cf5dd3
md""" 

У каждого диэлектрического материала существует так называемая длина волны Кристиансена, это точка, длина волны, на которой поверхность практически не отражает излучение ( следовательно в этой точке излучательная способность равна нулю) : 

Для выбранного материала **$( emissivity_type )** длина волна Кристиансена: 

``\lambda_{Ch}`` = **$( λ2[argmax(ϵ)])** μm

А излучательная способность в этой точке:

``\epsilon`` = **$(ϵ[argmax(ϵ)])** 

В соотвествии с формулой выше, при равенстве единице излучательной способности, пирометр будет регистрировать исключительно излучение поверхности, без паразитного отраженного излучения, так как коэффициент отражения будет равен нулю. 

"""

# ╔═╡ 5b647454-e5fc-4c65-8a0a-8eb499955cc1
if  use_custom
	p_custom = Pyrometers.Pyrometer( [custom_waves...] , type = :C  , ϵ = 1.0)
end

# ╔═╡ 12b14bc0-b3bb-48aa-b3fa-674b33b22034
@bind clear_plots Button("Clear plots")

# ╔═╡ aa7ddfa0-8e68-43c9-9051-30e31cd3d906
begin 
	clear_plots
	last_column_plot_data = Dict{String, Matrix{Float64}}()
end

# ╔═╡ 0428e504-e163-4945-a2fa-7175e0b33020
last_column_plot_data

# ╔═╡ ed653f67-2f49-432b-b691-a79f7e1d51e4
i_test = IsothermalSpectralQuantity(Base.Fix2(lamp_interpolator,240))

# ╔═╡ c9634225-fa52-4747-8f45-3511141bd164
md""" 
**Try plotly! : $(@bind is_use_plotly CheckBox(false))**

"""

# ╔═╡ efcc25ee-a507-4334-b338-f3e050a50c6b
use_custom && p_custom

# ╔═╡ e480137d-b6d9-4e18-92f0-640292bbb5f0
function two_planck(l , ϵ , T1 , T2)
	return ϵ * Planck.ibb.(l , T1) + (1.0 - ϵ) * Planck.ibb.(l , T2)
end

# ╔═╡ 302a58e0-ffba-464b-a962-7f85bca897c6
function generic_emission_spectrum(l , ϵ , T1 , U , lamp_interpolator)
	return ϵ * Planck.ibb.(l , T1) + (1.0 - ϵ) * lamp_interpolator(l , U)
end

# ╔═╡ cdf36fa2-9719-405c-86b9-2c78ee815dfe
function generic_emission_spectrum(l , ϵ , T1 , U , ::Nothing)
	return nothing
end

# ╔═╡ c1d897c7-377f-4a24-8694-963c4149b40b
function generic_emission_spectrum(l::AbstractVector , ϵ , T1 , U::Number, lamp_interpolator)
	return ϵ .* Planck.ibb.(l , T1) + (1.0 .- ϵ) .* lamp_interpolator.(l , U)
end

# ╔═╡ 6342e92b-4434-4e4b-aa2f-56405277caed
begin 
	I1 = @. ϵ * Planck.ibb(λ2 , T1)
	
	I2 = if !use_external_lamp 
		@. (1.0 - ϵ) * Planck.ibb.(λ2 , T2)
	else
		(1.0 .- ϵ) .* lamp_interpolator.(λ2 , T2)
	end
	I_measured =if !use_external_lamp 
		@. two_planck.(λ2 , ϵ , T1 , T2)
	else
		generic_emission_spectrum(λ2 , ϵ , T1 , T2 , lamp_interpolator)
	end
end;

# ╔═╡ 0404bf20-57a4-4c7c-bf23-70d3541a6787
begin 

	Tmeas = Dict{String , NamedTuple{(:Tmeas , :ΔT , :ΔTrel) ,Tuple{Float64 , Float64 , Float64}}}()
	if !isempty(selected_standard_pyros)
	for nm in selected_standard_pyros

		p = Pyrometers.Pyrometer(nm)
	
		is_two_wavelength_pyrometer = Pyrometers.is_spectral_band(p)
		# 
		l_cur = is_two_wavelength_pyrometer ? p.λ : [p.λ[]-0.2,p.λ[]+0.2 ]
		λ_pyr_interp = is_two_wavelength_pyrometer ? collect(range(l_cur...,length=30)) : p.λ

		e_evaluated = eint.(λ_pyr_interp)
		e_averaged = sum(e_evaluated)/length(e_evaluated)

		# setting averaged value of emissivity to pyrometer
		Pyrometers.set_emissivity!(p, e_averaged)
		
		# calculating the measured by the pyrometer value 
		I_cur = if !use_external_lamp
				two_planck.(λ_pyr_interp , e_evaluated , T1 , T2)
		else
			generic_emission_spectrum(λ_pyr_interp , e_evaluated , T1 , T2 , lamp_interpolator)
		end
		
		measure_intensity =  is_two_wavelength_pyrometer ? NumericalIntegration.integrate(λ_pyr_interp , I_cur) : I_cur[]
		# setting emissivity to one
		_tmeasured = Pyrometers.measure(p , measure_intensity, T_starting = T1)
		push!(Tmeas, String(p.type) =>(_tmeasured , _tmeasured - T1 , 100*(_tmeasured - T1)/T1))
	# temperature measured by the current pyrometer 
	end
	
	end
end

# ╔═╡ 01752ec9-b480-45b9-b299-d9c18d02a749
begin 
pyr_number = length(Tmeas)
pyr_data = Matrix{Any}(undef , (pyr_number,4))	
best_pyr = nothing
dT_opt = 1e6
for (i , (k , d)) in enumerate(Tmeas)
	pyr_data[i , 1] = k 
	pyr_data[i , 2] = d.Tmeas
	pyr_data[i , 3] = d.ΔT
	pyr_data[i , 4] = d.ΔTrel
	if d.ΔT < dT_opt
		global dT_opt = d.ΔT
		global best_pyr = Pair(k , d.Tmeas)
	end
end
pretty_table(HTML,pyr_data, column_labels= ["Тип","Тизм","ΔT, K","ΔT,%"] , title  ="Ошибка определения температуры поверхности пирометрами различных типов при температуре  поверхности образца T₁= $(T1) и U = $(T2) лампы $(interpolator_file_name)" )

end

# ╔═╡ 459f54a1-bbf0-4268-8bec-8142d436976a
begin 
	common_kwargs = (; fillrange=1e-1, fillalpha=0.3,dpi=600,legend_background_color=:white, legend_foreground_color = :black, legend_position=:right , grid = true, gridlinewidth=3, gridstyle = :dot,minorgrid=false, box = :on, linewidth = 3 , fmt = :png)
	
	(xmin , xmax) = extrema(λ_show)
	fl =@. (xmin <= λ2) & (λ2 <= xmax)
	y_lims = extrema(I_measured[fl])
	
	
	em_plot = Plots.plot(λ2 , eint.(λ2)  , label = nothing , title = emissivity_type; common_kwargs...)

	
	title_str = (emissivity_type == "fixed") ? " ϵ = $(ϵ)" :  "материала $(to_names[emissivity_type])"
	
	ppp = Plots.plot(λ2 , I_measured  , label = "Суммарное" ; common_kwargs...)
	Plots.plot!(ppp, λ2 , I1  , label = "Образец Tₛₐₘ. = $(T1)"  ; common_kwargs...)
	Plots.plot!(ppp, λ2 , I2  , label = "Лампа U = $(T2)" ; common_kwargs... , xscale=:log10 , yscale = :log10 , legend=:topright )
	title!(ppp , """Излучение $(title_str) 
		   Tbest = $(last(best_pyr)) ($(first(best_pyr)) - type) """)
	
	ylabel!(ppp , " Спектральная интенсивность , , W/m²⋅sr⋅μm")
	ylabel!(em_plot , " ϵ")
	ylims!(ppp , y_lims)
	ylims!(em_plot , (0.0,1.0))
	for p in (ppp , em_plot)
		xlabel!(p , " Длина волны , μm ")
		xlims!(p , (xmin , xmax))
		
	end
end

# ╔═╡ 5dafdc88-bd40-4347-aa62-e841d15c1bd7
em_plot

# ╔═╡ 18daa932-fd3a-4056-aa07-4dcf26c7d57a
ppp

# ╔═╡ 391ebf67-2bd9-4c28-a00d-7f3d732625c4
ppp

# ╔═╡ ac2343b5-6ea5-47b1-9c39-643cdf6d93af
begin 
	kvect = [k for k in keys(Tmeas)]
	(dT, ind_min) = findmin(k->Tmeas[k].ΔT , kvect)
	pyr_smaller_type = kvect[ind_min]
end;

# ╔═╡ 48b184a0-b1df-461e-a3f5-3c8be72ab875
md" **Выбор типа пирометра** : $(@bind selected_type Select(kvect , default = pyr_smaller_type) )"

# ╔═╡ 45396db0-b967-456a-b0e2-eec411675827
if is_recalculate #data preparation block

	T1_scan = 300.15:50:1773.15
	if !use_external_lamp
		T2_scan = 1273.15:50:3800.15
	else
		T2_scan = voltage_range
	end

	if use_custom_pyros
		p_selected = custom_selection
	elseif use_custom
		p_selected = p_custom
	else
		p_selected = filter(p->p.type == selected_type , pyrometers_vector2)[]
	end
	
	_is2 = Pyrometers.is_spectral_band(p_selected)
		# 
	l = _is2 ? p_selected.λ : [p_selected.λ[]-0.2 , p_selected.λ[]+0.2 ]
	
	_λ_pyr_interp = _is2 ? collect(range(l...,length=30)) : p_selected.λ

	e_s = eint.(_λ_pyr_interp)
	e_avg = sum(e_s)/length(e_s)

	# setting averaged value of emissivity to pyrometer
	Pyrometers.set_emissivity!(p_selected, e_avg)
	
	# calculating the measured by the pyrometer value 
	
	_N = length(T1_scan)
	_M = length(T2_scan)


end

# ╔═╡ 52312616-c471-4b35-bb98-46c42a82d191
begin 
	last_column_plot = nothing
	emissivity_type
	if !isempty(last_column_plot_data)
		last_column_plot = Plots.plot()
		for (k,d) in last_column_plot_data
			(k=="fixed") && continue
			Plots.plot!(last_column_plot , d[:,1] .- 273.15, d[:,2];plot_common_args..., label = to_names[k] , xticks =10  , linestyle = :auto , marker = :auto , markersize = 8 , markeralpha = 0.8)
		end
	end
	xlabel!(last_column_plot , "Температура поверхности, ᵒC")
	ylabel!(last_column_plot , "ΔT/T, %")
	title!(last_column_plot , Pyrometers.shorthand(p_selected))
end

# ╔═╡ 831ead94-fd87-4412-80e5-01f557eee355
begin 
	T_true = 957.98
	bb = Pyrometers.PlanckEmitter()
	e_surf = IsothermalSpectralQuantity(eint)
	refl = Pyrometers.SpectralReflectivity(e_surf)
	
	for u in (120.0 , 160.0 , 180.0 , 200.0 , 240.0)
		i_test = IsothermalSpectralQuantity(Base.Fix2(lamp_interpolator,u))
		i_full = e_surf * bb + refl * i_test
		T_source = p_selected(Pyrometers.fix_temperature( i_test , 1200.0) , e_surf)
		@show T_source
		i_fixed = Pyrometers.fix_temperature(i_full , T_true)
		Tmeasured = p_selected(i_fixed , e_surf)
		T_corrected = stray_radiation_corrected_temperature(p_selected ,Tmeasured, e_surf , T_source )
		@show T_corrected
		#T_corrected2 = stray_radiation_corrected_temperature(p_selected ,Tmeasured, e_surf , i_test)
	end
end

# ╔═╡ f1e1935e-5e1d-4cb7-a0b1-673989871786
Pyrometers.Planck.band_power(1200 , λₗ=p_selected.λ[1] , λᵣ=p_selected.λ[2])

# ╔═╡ 34502312-07df-4ce8-9bd5-0255c5741c4b
i_measured_check = Pyrometers.integrate(p_selected , 1200 , i_test)

# ╔═╡ 507e9cfa-4e66-47d6-83ed-d1d1e1de0676
T_corrected2 = stray_radiation_corrected_temperature(p_selected ,1234.5, e_surf , i_test)

# ╔═╡ 98b830a0-6505-4292-85ed-7c7aaa882bd9
Pyrometers.fit_integral(e_surf * bb , p_selected.λ[1] , p_selected.λ[2] , i_measured_check)

# ╔═╡ 3e19d251-91f6-4383-bfc8-ffe816570f42
if is_recalculate
	md"""
	Tsurface = $(@bind T_surf_selected Select(T1_scan , default = 1273.15))
	
	Tlamp = $(@bind T_lamp_selected Select(T2_scan , default = 2873.15))
	
	"""
end	

# ╔═╡ 10298d52-d411-475f-b7f6-8562ed2a25bc
if is_recalculate 
	

	
	ΔT_mat = zeros((_M , _N))
	Tmeas_mat = zeros((_M , _N))
	Pyrometers.set_emissivity!(p_selected , 1.0)
	Threads.@sync for j in 1 : _N
		Threads.@spawn begin
		for i in 1 : _M
			t1 =  T1_scan[j]
			t2 = T2_scan[i]
			
			i_measured_selected = if !use_external_lamp
				two_planck.(_λ_pyr_interp , e_s , t1, t2)
			else
				generic_emission_spectrum(_λ_pyr_interp , e_s , t1 , t2 , lamp_interpolator)
			end
			
			i_measured_int =  _is2 ? NumericalIntegration.integrate(_λ_pyr_interp , i_measured_selected) : i_measured_selected[]

			
			
			t_meas = Pyrometers.measure(p_selected , i_measured_int, T_starting = t1)
			
			ΔT_mat[i , j] = (t_meas - t1)/t1
			Tmeas_mat[i , j] = t_meas
		end
		end
	end
	
end

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
		#=p_res = Plots.surface(T1_scan .-273.15, T2_scan , 100.0*ΔT_mat )
		xlabel!("Температура поверхности, ᵒС")
		ylabel!("Напряжение на лампе, В")
		zlabel!("ΔT/T, %")=#
		 p_res = Plots.surface(
        T1_scan .- 273.15, 
        T2_scan, 
        100.0 .* ΔT_mat,
        
        # Цветовая палитра и сглаживание
        colormap = :viridis,           # Аналог цвета из вашего скрипта Plotly
        colorbar = true,               # Отображение шкалы цветов
        linealpha = 0.5,               # Делает линии сетки на поверхности тонкими/полупрозрачными
        fillalpha = 0.95,              # Небольшая прозрачность для мягкости теней
        
        # Подписи осей и шрифты
        xlabel = "Температура поверхности, ᵒС \n",
        ylabel = "   Напряжение на лампе, В\n",
        zlabel = "ΔT/T, %",
        guidefontsize = 14,            # Размер шрифта подписей осей
        tickfontsize = 10,              
        
        # Размер и пропорции
        size = (800, 600),
        camera = (45, 30), 
		title = "$(Pyrometers.shorthand(p_selected)) : $(to_names[emissivity_type])",
			 fmt = :png
    )
	end
end

# ╔═╡ ff5266ce-b7f2-46e2-a4b6-30057879bd97
last_column_plot_data[emissivity_type] = hcat(collect(T1_scan) , abs.(vec(100*ΔT_mat[end , :])));

# ╔═╡ ce4f2fdd-16b1-46e8-88a4-a952896b6df8
if is_recalculate

	j_selected = findfirst(t->t == T_surf_selected , T1_scan)
	i_selected = findfirst(t->t == T_lamp_selected , T2_scan)
	md" Tmeas ± ΔT = $( Tmeas_mat[i_selected , j_selected] ) ± $(ΔT_mat[i_selected , j_selected] *  T1_scan[j_selected] )"
end

# ╔═╡ dd1561e2-233f-425a-832f-130b49f0bf0b
if is_recalculate
	
	out_table =map(t->round(t,digits=2),  hcat(T2_scan, 100*ΔT_mat[: , 1:4:end]))
	col_nms =vcat("V",map(t -> "$(round(t - 273.15,digits =1))",  T1_scan[1:4:end]))
	pretty_table(HTML , out_table , title="$(to_names[emissivity_type])", column_labels=col_nms)
end

# ╔═╡ 17131508-d23b-4e95-9303-e4b685dda736
maximum(ΔT_mat)

# ╔═╡ 05ec0ea2-cfec-4e91-a46a-68bbdaefd562
if is_recalculate
	out_table_T = vcat(hcat([0.0] , transpose(T1_scan)) , hcat(T2_scan, Tmeas_mat))
	pretty_table(HTML , out_table_T )
end

# ╔═╡ Cell order:
# ╠═30743a02-c643-4bdc-837e-b97299f9520a
# ╠═5e712312-0fc7-4205-84cc-834d57b814a3
# ╠═ba23c985-74c4-41f3-8bc4-f7287e30e47f
# ╠═8911823a-2210-4328-b9a5-aa88af2146c5
# ╠═b1e73e8f-8eb4-4688-86a8-b787d6b15454
# ╠═89a11dcd-b3b5-4349-930d-a66ad74e8fa2
# ╠═4f1397e1-018d-4965-b3e1-78ccac62b54c
# ╠═171409eb-22b5-4bc5-a8e2-eac0932a24f3
# ╠═04f9f75e-34cc-4921-971e-6958ed960518
# ╠═7d37ae06-c81b-4aa1-a05c-e03b83e395d7
# ╟─7cfee16a-4c47-4ded-9b5a-4a21e59dd9ad
# ╟─d5ee3913-66be-47d7-a755-699ba64b4f98
# ╟─d442014a-20e6-4be4-ac7f-f13de329dec5
# ╟─27b3c586-9eb0-4a51-b9ca-a9c0379fccdf
# ╟─f22d22b6-5d98-4cc4-998f-a53e92809618
# ╟─51de94a8-60ae-4099-97f8-f20582c4bc48
# ╟─72095f49-f6c0-41a8-be2f-323a0bacfc51
# ╟─7cc110e9-7655-4dfc-b1e0-ab3905866425
# ╟─4d6337aa-cfc7-4154-a395-5aa53e23d01a
# ╟─f5a2937c-ce1c-4c03-8e86-8db363bd21c9
# ╟─e64ba0f6-48a6-4e87-a347-e824f350ef4a
# ╠═1e429893-8aaa-4e95-aff8-4eab826a5c4f
# ╟─03d76e64-ebf4-432b-b9be-d4cb26275f55
# ╟─2ad3ec82-54a2-49ac-94ef-579f808dfb1a
# ╟─8a066ee5-80e9-462f-9a61-15851468aa63
# ╟─b7fac177-c211-4635-992f-e6473be7bdae
# ╟─69f9052b-d925-49fb-8eb9-089671e64339
# ╟─e2a9aa39-2490-4681-89d3-a01f058f6feb
# ╟─02eee968-ff43-4b82-8d68-efede1a220dd
# ╟─6c2d1616-4992-4123-a356-0545cf81ce08
# ╟─9b08b767-7e8f-4483-9f2f-226022ce10e4
# ╟─c69acbf6-94fb-4ac3-8d56-d1f9dda11440
# ╟─7071a6f4-e296-4e53-8e6f-24f1f038c1a5
# ╟─712828a7-fb54-42e6-95fc-233243190f59
# ╟─f763d449-2a7a-4008-a183-823a774bc25e
# ╟─667f7c30-56e0-461f-b35b-c924007eb9f2
# ╟─a861d56f-f6c9-4754-b7a9-ed63713f1f2f
# ╟─d08ec8f2-7689-4043-9f28-da06ab0124b9
# ╟─3ddc9d4e-ff96-4173-899b-3eae0facb830
# ╟─c5ac80ee-8143-4c28-bbff-2ac761c71fac
# ╟─0c9fe7b1-374c-4fb8-9cfe-9337389713bf
# ╟─bc2d93ae-6c30-462c-96e0-30fdb84d7c63
# ╟─d3199b6e-9779-4def-b701-fe85d1035045
# ╟─6aa82819-6f5a-43af-ac03-4775fce51564
# ╟─72947d97-0a97-4064-a2fe-08d19dec0f0e
# ╟─854f8c06-aac6-464f-899c-41ca706b259a
# ╟─9ce196c1-8915-46da-9aba-f13d7655959a
# ╟─86af6afc-b28a-4e84-952a-bd29710374f8
# ╟─efc35420-d0e6-4795-94b6-d43289b4de44
# ╟─6342e92b-4434-4e4b-aa2f-56405277caed
# ╟─620dce38-97ce-495f-9b23-1b8290cbd973
# ╟─b15bd3b0-0971-4b22-a109-0ed56afda615
# ╠═da0a4a00-692a-4d29-97cb-74366e8374c9
# ╟─90e15591-d885-48f3-bc46-cb9385c9ed13
# ╟─afb61470-2455-4742-b81f-63baf5a4c3ce
# ╟─8b441a27-004a-4468-b89f-f3a6cc6f9bbf
# ╠═0f0baa19-aa64-4daa-a501-861e22cdafde
# ╠═f3373768-d598-4d1e-b01b-1929e233f36a
# ╠═6edddd01-915c-45bb-b3ef-155f4dcd9d8c
# ╟─bfce276c-389c-405c-a552-979ed38133a2
# ╟─31a59906-e0d8-4ae8-b1ca-34b9f01062c8
# ╠═b27eff45-93e2-4aeb-a7ef-88132dd7e9e8
# ╠═0cb1cfa0-6f2d-49bc-ac35-bb76688bb510
# ╠═4cc48476-fb0c-4a30-a93d-2e6caa53e232
# ╠═ac15bdd0-7679-4511-b0be-4885f2d366bd
# ╠═d67ef942-be2d-4fad-91ba-eff9eda9d8e3
# ╟─5dafdc88-bd40-4347-aa62-e841d15c1bd7
# ╟─a7ff8a2d-a81d-4474-b27f-565de2cf5dd3
# ╟─18daa932-fd3a-4056-aa07-4dcf26c7d57a
# ╟─8cef05a1-2974-4c38-b73e-fa706f347fcc
# ╟─36ba2396-bb5e-4d22-a58e-9ab27cd18b2d
# ╟─91bbd553-4e4a-431d-9d54-b0f4882fd426
# ╟─15f1519b-d924-4fa9-b212-eebba75c544a
# ╟─1c2bd109-8023-49a5-87ba-32d755b201dd
# ╟─01752ec9-b480-45b9-b299-d9c18d02a749
# ╟─0404bf20-57a4-4c7c-bf23-70d3541a6787
# ╠═ad622962-9cd8-432c-a988-ec8b3d676077
# ╟─c918b685-1b6e-4ce3-a1e6-501c6b4ffe0e
# ╠═96feae26-99e4-45c9-869b-410ec7965b63
# ╠═c5a00df1-ff02-4564-9838-e678f8ea7389
# ╟─459f54a1-bbf0-4268-8bec-8142d436976a
# ╠═391ebf67-2bd9-4c28-a00d-7f3d732625c4
# ╟─1811e43c-f7db-47b1-9b83-bb38455d7db3
# ╟─ac2343b5-6ea5-47b1-9c39-643cdf6d93af
# ╟─8ece7476-4257-49f8-a400-50bcc347908d
# ╟─54339700-71fd-48bf-a2ef-0c3267b9d81b
# ╟─7f76bc22-77c3-4eb3-9d49-588e653df2e7
# ╟─a0197a9a-34bf-4a3e-af8a-c23ea777f482
# ╟─a90bde8a-222b-4b93-892b-9fa4cfc43372
# ╟─d86d7e8d-94e7-4d1b-b4d8-df52e338935d
# ╟─7793976e-6714-4bc1-9d97-412dd2a67480
# ╟─c01e05c7-2b29-4f57-96ee-9a4eebef9025
# ╟─48b184a0-b1df-461e-a3f5-3c8be72ab875
# ╟─f181980f-bf72-4468-8daa-9461c6c901e0
# ╠═0cb50b02-ab41-416c-8610-c3ff318b117b
# ╟─5b647454-e5fc-4c65-8a0a-8eb499955cc1
# ╟─45396db0-b967-456a-b0e2-eec411675827
# ╟─12b14bc0-b3bb-48aa-b3fa-674b33b22034
# ╟─aa7ddfa0-8e68-43c9-9051-30e31cd3d906
# ╟─ff5266ce-b7f2-46e2-a4b6-30057879bd97
# ╟─0428e504-e163-4945-a2fa-7175e0b33020
# ╟─52312616-c471-4b35-bb98-46c42a82d191
# ╠═10298d52-d411-475f-b7f6-8562ed2a25bc
# ╠═831ead94-fd87-4412-80e5-01f557eee355
# ╠═f1e1935e-5e1d-4cb7-a0b1-673989871786
# ╠═ed653f67-2f49-432b-b691-a79f7e1d51e4
# ╠═34502312-07df-4ce8-9bd5-0255c5741c4b
# ╠═507e9cfa-4e66-47d6-83ed-d1d1e1de0676
# ╠═98b830a0-6505-4292-85ed-7c7aaa882bd9
# ╟─c9634225-fa52-4747-8f45-3511141bd164
# ╟─3e19d251-91f6-4383-bfc8-ffe816570f42
# ╟─ce4f2fdd-16b1-46e8-88a4-a952896b6df8
# ╟─efcc25ee-a507-4334-b338-f3e050a50c6b
# ╟─dd1561e2-233f-425a-832f-130b49f0bf0b
# ╠═17131508-d23b-4e95-9303-e4b685dda736
# ╟─05ec0ea2-cfec-4e91-a46a-68bbdaefd562
# ╟─e480137d-b6d9-4e18-92f0-640292bbb5f0
# ╟─302a58e0-ffba-464b-a962-7f85bca897c6
# ╟─cdf36fa2-9719-405c-86b9-2c78ee815dfe
# ╠═c1d897c7-377f-4a24-8694-963c4149b40b
