const year = 365days

## Surface PAR and turbulent vertical diffusivity based on idealised mixed layer depth 
@inline PAR⁰(x, y, t) =
    60 *
    (1 - cos((t + 15days) * 2π / year)) *
    (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2
