[![Build Status](https://github.com/Manarom/Pyrometers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Manarom/Pyrometers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://manarom.github.io/Pyrometers.jl)



# Pyrometers.jl

`Pyrometers.jl` is a Julia package designed for modeling, simulating, and evaluating non-contact temperature measurements focused on preformance. It provides calculation methods for brightness (radiation) and two-color (ratio) pyrometers operating at discrete wavelengths or over integrated spectral bands.

The package is engineered for zero-allocation execution, making it optimal for real-time control loops, industrial data parsing, and automatic differentiation.

---

##  Main Functionality

* **Temperature Inversion (`measure`):** Solves for the actual surface temperature from a given raw detector signal power handle using fast root-finding methods.
* **Reflection Compensation (`corrected_temperature`):** Subtracts ambient background radiation, stray furnace reflections, and grey-body noise from raw measurements to find the reflection-free target temperature.
* **Signal Generation (`signal`):** Simulates the expected radiant power or intensity value hitting a detector for any given target temperature.
* **Calibration Optimization (`fit_ϵ` / `fit_ϵ!`):** Back-calculates the necessary effective emissivity or spectral ratio slope (`e_slope`) required to match a pyrometer's reading with a known reference temperature.

There are several default pyrometers with specified spectral ranges, see the following figure:

<p float="left">
  <img src="./notebooks/Pyrometers.png" width="400"/>
</p>

---

## 🚀 Quick Start

### 1. The Core Temperature Solver
Feed an incoming raw detector signal power to your pyrometer, and get the actual temperature back instantly:

```julia
using Pyrometers 
import Pyrometers.Planck as PF # package for thermal radiation is PlanckFunction.jl
ϵ = 0.55 # surface emissivity
T_real = 1234.89
# Create a spectral band pyrometer (2.4 μm to 4.5 μm) calibrated for ϵ = 0.55
p = SpectralBandPyrometer(2.4, 4.5, ϵ = ϵ ) 
# Create a two-color  band pyrometer (2.4 μm to 4.5 μm) and (6.7 μm to 9.0 μm) 
p_color = TwoBandsRatioPyrometer((2.4 , 4.5) , (6.7 , 9.0) , ϵ1 = ϵ , ϵ2 = 0.2 ) 
i = PF.band_power(T_real ,  λₗ = 2.4 , λᵣ = 4.5) # total intensity in pyrometer spectral range
r = PF.spectral_band_ratio((2.4 , 4.5) , (6.7 , 9.0) , T_real) # evaluates to band power ratio
r *= ( ϵ /0.2)
# Calculate temperature using the shortcut functor syntax p(signal value)
T_measured = p( ϵ * i ) 
T_measured2 = p_color(r)
println("Treal : $(T_real) , K")
println("Band pyrometer: $(T_measured) K")
println("Ratio pyrometer: $(T_measured2) K")

```

### 2. De-Noising Ambient Furnace Reflections
If your target is inside a hot oven, background reflections distort your readings. Wipe them out with `corrected_temperature`:

```julia

  e_real = 0.35
  p_band = SpectralBandPyrometer(2.4, 2.5 ,  ϵ = e_real)
  T_real = 900.11 
  T_stray = 1000.0 # temperature of stray radiation
  isurface = e_real * PF.band_power( T_real,  λₗ = 2.4 , λᵣ = 2.5) # real temperature is 1700.11
  reflected = (1 - e_real) * PF.band_power( T_stray,  λₗ = 2.4 , λᵣ = 2.5)
  i = isurface + reflected
  T_measured = p_band(i) # measured temperature is  incorrect due to stray radiation
  T_corrected = corrected_temperature(p_band , T_measured ,T_stray , 1.0) # last argument is surrounding effective emissivity
  # Corrected temperature is equal to T_real
```

### 3. Calculating Integral Emissivity
Easily average spectral emissivity over the Planck distribution using either custom formulas or raw data tables:

```julia
p_band = SpectralBandPyrometer(2.4, 4.5)

# Approach A: Using an analytical formula or spline function
ϵ_formula(λ) = 0.45 + 0.02 * λ
ϵ_avg = integral_emissivity(p_band, ϵ_formula, 1273.15)

# Approach B: Using raw experimental data tables
λ_grid = range(1.0 , 10.0 , 100)
ϵ_grid = ϵ_formula.(λ_grid)

(ϵ_avg,)= integral_emissivity(p_band, λ_grid, ϵ_grid, 1273.15)
```
### 4. Converting temperatures 
Estimating the emissivity of the surface from the temperature measured by the pyrometer and 
the real temperature 

```julia

p_band = SpectralBandPyrometer(2.4, 4.5 ,  ϵ=1.0)
e_real = 0.35
T_real = 1700.11
i = e_real * PF.band_power( T_real,  λₗ = 2.4 , λᵣ = 4.5) # real temperature is 1700.11
T_measured = p_band(i) # measured temperature is for incorrect emissivity
T_corrected = convert_temperature(p_band , T_measured , e_real) # real emissivity is 0.35
# Corrected temperature is equal to T_real
fit_ϵ!(p_band , T_measured , T_real) # fits the emissivity of the pyrometer 

```

---

## 📄 License

This package is open-source software licensed under the [MIT License](LICENSE).

