[![Build Status](https://github.com/Manarom/Pyrometers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Manarom/Pyrometers.jl/actions/workflows/CI.yml?query=branch%3Amain)

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://manarom.github.io/Pyrometers.jl)
# Pyrometers.jl

`Pyrometers.jl` is a Julia package designed for modeling, simulating, and evaluating non-contact temperature measurements. It provides fast calculation methods for brightness (radiation) and two-color (ratio) pyrometers operating at single wavelengths or over integrated spectral bands. It also supports multiwavelength pyrometers operating at multiple wavelengths simultaneously, allowing the decoupling of emissivity from blackbody wavelength dependence using non-linear least-squares fitting.

The package is engineered from the ground up for strict **zero-allocation execution** and full type stability.

The purpose of this package is to simulate virtual pyrometers of **five different types**:

* `SpectralBandPyrometer`: A classical pyrometer operating within a fixed wide spectral band.
* `TwoBandsRatioPyrometer`: A spectral ratio pyrometer that accounts for the finite width of two channels.
* `SingleWavelengthPyrometer`: A standard pyrometer operating at a single, fixed wavelength.
* `TwoWavelengthRatioPyrometer`: A classical spectral ratio pyrometer that does not account for the finite width of its spectral channels.
* `MultiwavelengthPyrometer`: A pyrometer based on least-squares fitting of discrete spectral intensity data across multiple wavelengths (theoretically eliminating the need to know surface emissivity).

---

## Main Functionality

* **Temperature Inversion (`measure` / Functor Syntax):** Solves for the actual surface temperature from a given raw detector signal power or discrete intensity using high-order root-finding methods (Halley's method) or least square fitting using Levenberg-Marquardt method (custom implementation on StaticArrays) for multiwavelength pyrometer.
* **Stray Radiation Compensation (`stray_radiation_corrected_temperature` and `external_source_corrected_temperature`):** Isolates and removes parasitic reflected background radiation, hot furnace wall reflections, and gray-body cavity noise from raw measurements to extract the true target temperature.
* **Signal Generation (`signal`):** Simulates  radiant power or intensity value hitting a detector for any given target temperature and emissivity configuration.
* **Lazy Spectral Algebra:** Features an embedded symbolic-numeric composition engine (`*`, `/`, `+`) for mixing complex temperature-dependent spectral functions (`AbstractSpectralQuantity`) and raw tabular data matrices (`AbstractDiscreteQuantity`) with zero allocation overhead.
* **Calibration Optimization (`fit_ϵ` / `fit_ϵ!`):** Back-calculates the effective emissivity or spectral ratio slope required to match a pyrometer's reading with a known reference temperature.

Several pre-defined industrial pyrometers with specific spectral ranges are available out of the box (see the figure below):

<p float="left">
  <img src="./notebooks/Pyrometers.png" width="400" alt="Pyrometers Spectral Ranges"/>
</p>

---

## Quick Start

### 1. Core Temperature Solver
Feed an incoming raw detector signal power or ratio handle to your pyrometer, and get the actual temperature back instantly:

```julia
using Pyrometers 
import Pyrometers.Planck as PF # The package for thermal radiation is PlanckFunction.jl

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

# Multiwavelength pyrometer
l = range(1, 2, 30)
i = 0.3 * PF.ibb.(l, 1273.15) # thermal emission spectrum of a surface with constant emissivity 
p_multiwavelength = MultiwavelengthPyrometer{30}(l) # by default approximates emissivity with a Bernstein polynomial
T = p_multiwavelength(i) # evaluates temperature without a priori emissivity knowledge
# T ≈ 1273.15
e = emissivity_poly(p_multiwavelength) # returns fitted emissivity polynomial approximation as a callable scaled polynomial object
e.(1:0.1:1.5) # returns emissivity evaluated at new wavelengths
```

### 2. De-Noising Ambient Furnace Reflections (Advanced Radiosity)
If your target is inside a hot oven or narrow cavity, background reflections distort your readings. Wipe them out by evaluating mutual reflection configurations (such as an enclosure cavity or finite parallel plates) using geometric factors (ξ):

```julia
using Pyrometers 

ϵ_surf = IsothermalSpectralQuantity(l -> 0.6 + l/20) # surface emissivity 
ϵ_wall = IsothermalSpectralQuantity(l -> 0.9 - l/20) # external source emissivity 
Tsource = 1500.0 # source temperature
Ttrue = 987.5 # true temperature of the surface

# Imitating an external radiation source 
bb = PlanckEmitter() # creating external source imitator 
i_incident = fix_temperature(bb, Tsource)
refl = SpectralReflectivity(ϵ_surf)

p = SpectralBandPyrometer(2.0, 3.0)
i_full = ϵ_surf * bb + refl * ϵ_wall * i_incident # spectral algebra usage
i_full_iso = fix_temperature(i_full, Ttrue)
I_total = integrate(p, i_full_iso)
T_meas = p(I_total, ϵ_surf) # measured temperature including stray radiation impact

geom = EnclosureGeometry()
T_corrected = external_source_corrected_temperature(p, T_meas, ϵ_surf, Tsource, ϵ_wall, geom) # applying correction: T_corrected ≈ Ttrue
```

### 3. Combining Lazy Spectral Quantities with Algebra
Build complex custom multi-layer or selective emission models on the fly using native algebraic operators. The underlying code evaluates the necessary first- and second-order derivatives (`eval_Dₜ`) over the composite chain rule without a single heap allocation:

```julia
import Pyrometers as P

i_source = P.PlanckEmitter() # plack spectrum emitter 
e_source = P.IsothermalSpectralQuantity(l -> 0.9 + 1e-2 * l)

# Create a custom temperature-dependent compound profile using lazy products
i_corrected_source = i_source * e_source 

p = P.SpectralBandPyrometer(2.0, 3.0)
P.set_emissivity!(p, 0.6)

# Evaluates the full model composition and filters background noise safely
T_clean = P.stray_radiation_corrected_temperature(p, 1200.0, i_corrected_source, 1300.0)

# Multiwavelength pyrometry measurement 
p_multi = P.MultiwavelengthPyrometer{50}(range(1, 2, 50))
i = P.fix_temperature(i_source * e_source, 1345.6) # fixing source true temperature 
T = p_multi(i) # T ≈ 1345.6
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

## Roadmap & Future TODOs

The long-term objective of this ecosystem is to establish a comprehensive, unified pyrometry toolkit in Julia. The immediate development focus is centered on merging classical radiative instruments with data-driven spectral inversion techniques.

* Add stray radiation and external source exclusion for multiwavelength pyrometers.
* Add geometrical integration of external radiation (through multiple view factors).

---

## License

This package is open-source software licensed under the [MIT License](LICENSE).