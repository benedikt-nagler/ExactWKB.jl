```@meta
CurrentModule = ExactWKB
```

# Plotting

Plotting lives in a package extension, loaded on demand. Calling one of these functions before
a Makie backend is loaded raises an error telling you what to load.

```julia
using CairoMakie                          # or GLMakie

g = stokes_graph(prob; theta = 0.3)
plot_stokes_graph(g)                      # turning points, traced lines, saddles

plot_stokes_graph(stokes_graph_family(prob, [-0.15, 0.0, 0.15]))   # a θ-family, side by side

t = ideal_triangulation(g)
plot_triangulation(g, t)                  # the dual triangulation over the graph
```

The Borel-plane plots of the Voros series come from `Resurgence.jl` and work on anything this
package produces:

```julia
import Resurgence
Resurgence.plot_borel_plane(
    Resurgence.pade(Resurgence.borel(quantum_series(vs)); reduce = true))
```

```@docs
plot_stokes_graph
plot_triangulation
```
