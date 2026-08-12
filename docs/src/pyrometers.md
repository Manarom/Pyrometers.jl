# Pyrometers.jl

`Pyrometers.jl` is a Julia package designed for modeling, simulating, and evaluating non-contact temperature measurements focused on preformance. It provides calculation methods for brightness (radiation) and two-color (ratio) pyrometers operating at discrete wavelengths or over integrated spectral bands.

The package is engineered for zero-allocation execution, making it optimal for real-time control loops, industrial data parsing, and automatic differentiation.


```@autodocs
    Modules = [Pyrometers]
    Order   = [:module,:type,:constant,:function]
```