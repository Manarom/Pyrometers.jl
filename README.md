[![Build Status](https://github.com/Manarom/Pyrometers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Manarom/Pyrometers.jl/actions/workflows/CI.yml?query=branch%3Amain)

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://manarom.github.io/Pyrometers.jl)

# Pyrometers.jl

`Pyrometers.jl` is a high-performance Julia package designed for modeling, simulating, and evaluating non-contact temperature measurements. It provides fast calculation methods for brightness (radiation) and two-color (ratio) pyrometers operating at single wavelengths or over integrated spectral bands.

The package is engineered from the ground up for strict **zero-allocation execution** and full type stability, making it ideal for real-time control loops, high-frequency industrial data parsing, and automatic differentiation.

---

## 🛠 Main Functionality

* **Temperature Inversion (`measure` / Functor Syntax):** Solves for the actual surface temperature from a given raw detector signal power handle using high-order root-finding methods (Halley's method).
* **Stray Radiation Compensation (`stray_radiation_corrected_temperature`):** Isolates and removes parasitic reflected background radiation, hot furnace wall reflections, and grey-body cavity noise from raw measurements to isolate the true target temperature.
* **Signal Generation (`signal`):** Simulates the exact radiant power or intensity value hitting a detector for any given target temperature and emissivity configuration.
* **Lazy Spectral Algebra:** Features an embedded symbolic-numeric composition engine (`*`, `/`, `+`) for mixing complex temperature-dependent spectral functions (`AbstractSpectralQuantity`) and raw tabular data matrices (`AbstractDiscreteQuantity`) with zero allocation overhead.
* **Calibration Optimization (`fit_ϵ` / `fit_ϵ!`):** Back-calculates the necessary effective emissivity or spectral ratio slope required to match a pyrometer's reading with a known reference temperature.

There are several pre-defined industrial pyrometers with specific spectral ranges available out of the box (see the figure below):

<p float="left">
  <img src="./notebooks/Pyrometers.png" width="400"/>
</p>

---

## 🚀 Quick Start

### 1. Core Temperature Solver
Feed an incoming raw detector signal power or ratio handle to your pyrometer, and get the actual temperature back instantly:

```julia
using Pyrometers 
import Pyrometers.Planck as PF # package for thermal radiation is PlanckFunction.jl

ϵ = 0.55 # surface emissivity
T_real = 1234.89

# Create a spectral band pyrometer (2.4 μm to 4.5 μm) calibrated for ϵ = 0.55
p = SpectralBandPyrometer(2.4, 4.5, ϵ = ϵ) 

# Create a two-color band ratio pyrometer (2.4-4.5 μm) and (6.7-9.0 μm) 
p_color = TwoBandsRatioPyrometer((2.4, 4.5), (6.7, 9.0), ϵ1 = ϵ, ϵ2 = 0.2) 

i = PF.band_power(T_real, λₗ = 2.4, λᵣ = 4.5) # total intensity in pyrometer spectral range
r = PF.spectral_band_ratio((2.4, 4.5), (6.7, 9.0), T_real) # evaluates to band power ratio
r *= (ϵ / 0.2)

# Calculate temperature using the zero-allocation shortcut functor syntax p(signal value)
T_measured = p(ϵ * i) 
T_measured2 = p_color(r)

println("T_real: \$(T_real) K")
println("Band pyrometer: \$(T_measured) K")
println("Ratio pyrometer: \$(T_measured2) K")
```

### 2. De-Noising Ambient Furnace Reflections (Advanced Radiosity)
If your target is inside a hot oven or narrow cavity, background reflections distort your readings. Wipe them out by evaluating mutual reflection configurations (such as an enclosure cavity or finite parallel plates) using geometric factors (ξ):

```julia
using Pyrometers

p = SpectralBandPyrometer(2.0, 4.5)
T_measured = 1200.0  # Raw contaminated temperature reading
T_furnace  = 1500.0  # Temperature of the heated parasitic source/walls

# Define arbitrary temperature-dependent spectral profiles for your target and furnace walls
# Using the built-in automatic differentiation backend 
ϵ_object = GenericDifferentiableSpectralQuantity((λ, t) -> 0.6 - 0.0001 * t)
ϵ_wall   = 0.85 # Supports constant gray numbers seamlessly too

# Case A: Small target enclosed inside an enormous furnace cavity (ξ = 0.0)
T_true_enc = stray_radiation_corrected_temperature(p, T_measured, ϵ_object, T_furnace, ϵ_wall, EnclosureGeometry())

# Case B: Closely spaced comparable surfaces or infinite parallel plates (ξ = 1.0)
T_true_par = stray_radiation_corrected_temperature(p, T_measured, ϵ_object, T_furnace, ϵ_wall, ParallelGeometry())

# Case C: Custom geometry using an area-weighted view factor ξ = F₁₂ * A₁ / A₂
geom = ViewFactorGeometry(0.45)
T_true_custom = stray_radiation_corrected_temperature(p, T_measured, ϵ_object, T_furnace, ϵ_wall, geom)
```

### 3. Combining Lazy Spectral Quantities with Algebra
Build complex custom multi-layer or selective emission models on the fly using native algebraic operators. The underlying code evaluates the necessary first and second-order derivatives (`eval_Dₜ`) over the composite chain rule without a single heap allocation:

```julia
import Pyrometers as P

i_source = P.PlanckEmitter() 
e_source = P.IsothermalSpectralQuantity(l -> 0.9 + 1e-2 * l)

# Create a custom temperature-dependent compound profile using lazy products
i_corrected_source = i_source * e_source 

p = P.SpectralBandPyrometer(2.0, 3.0)
P.set_emissivity!(p, 0.6)

# Evaluates the full model composition and filters background noise safely
T_clean = P.stray_radiation_corrected_temperature(p, 1200.0, i_corrected_source, 1300.0)
```

### 4. Converting Operating Temperatures & Emissivity Tuning
Quickly transpose temperature scales between alternative material calibrations or align a device's baseline parameters to match a rigorous thermal validation standard:

```julia
p_band = SpectralBandPyrometer(2.4, 4.5, ϵ = 1.0)
e_real = 0.35
T_real = 1700.11

i = e_real * PF.band_power(T_real, λₗ = 2.4, λᵣ = 4.5)
T_measured = p_band(i) # Incorrect temperature due to blackbody scale assumption (ϵ = 1.0)

# Translate the reading back to a scale modeled after a true emissivity of 0.35
T_corrected = convert_temperature(p_band, T_measured, e_real) 

# Mutate and fit the baseline internal pyrometer calibration profile to match reference data
fit_ϵ!(p_band, T_measured, T_real) 
```

---


## 📌 Roadmap & Future TODOs

The long-term objective of this ecosystem is to establish a comprehensive, unified pyrometry toolkit in Julia. The immediate development focus is centered on merging classical radiative instruments with data-driven spectral inversion techniques.

* **Integrate Multiwavelength Capabilities:** Port and combine the core inversion algorithms from [BandPyrometry.jl](https://github.com/Manarom/BandPyrometry.jl.git) into this package.
* **Emissivity-Agnostic Pyrometry:** Incorporate advanced multiwavelength pyrometry methods capable of reconstructing target surface temperatures from continuous thermal emission spectra *without requiring an a priori known emissivity profile*.
* **Unified API Design:** Blend classical brightness/ratio equations with constrained polynomial optimization (e.g., via `ScaledPolynomials.jl`) into a single, cohesive interface.
* **Maintain Zero-Allocation Performance:** Ensure that incoming multiwavelength optimization routines and data-parsing pipelines respect the performance constraints of the package core (strict type-stability, static matrix sizing, and heap-allocation-free calculations).

---

## 📄 License

This package is open-source software licensed under the [MIT License](LICENSE).
