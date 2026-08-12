`Pyrometers.jl` is a Julia package designed for modeling, simulating, and evaluating non-contact temperature measurements focused on preformance. It provides calculation methods for brightness (radiation) and two-color (ratio) pyrometers operating at discrete wavelengths or over integrated spectral bands.

The package is engineered for zero-allocation execution, making it optimal for real-time control loops, industrial data parsing, and automatic differentiation.


# Pyrometers.jl

`Pyrometers.jl` is a Julia package designed for modeling, simulating, and evaluating non-contact temperature measurements focused on preformance. It provides calculation methods for brightness (radiation) and two-color (ratio) pyrometers operating at discrete wavelengths or over integrated spectral bands.

The package is engineered for zero-allocation execution, making it optimal for real-time control loops, industrial data parsing, and automatic differentiation.

---

##  Main Functionality

* **Temperature Inversion (`measure`):** Solves for the actual surface temperature from a given raw detector signal power handle using fast root-finding methods.
* **Reflection Compensation (`corrected_temperature`):** Subtracts ambient background radiation, stray furnace reflections, and grey-body noise from raw measurements to find the reflection-free target temperature.
* **Signal Generation (`signal`):** Simulates the expected radiant power or intensity value hitting a detector for any given target temperature.
* **Calibration Optimization (`fit_ϵ` / `fit_ϵ!`):** Back-calculates the necessary effective emissivity or spectral ratio slope (`e_slope`) required to match a pyrometer's reading with a known reference temperature.