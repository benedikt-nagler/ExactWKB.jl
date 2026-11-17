```@meta
CurrentModule = ExactWKB
```

# Stokes graphs

The Stokes graph at angle ``\theta`` is the level set

```math
\mathrm{Im}\left(e^{-i\theta}\int_{z_*}^{z}\sqrt{Q}\,dz\right) = 0
```

emanating from each turning point ``z_*``. It records which Borel-summed solution you are on:
crossing a Stokes line changes the dominant/recessive labelling, and as ``\theta`` rotates the
graph's *topology* changes only at finitely many angles.

[`stokes_graph`](@ref) traces the lines in an augmented state ``(z, w)`` with
``w = \sqrt{Q}`` carried as a dynamical variable,

```math
\frac{dz}{dt} = \frac{e^{i\theta}}{w}, \qquad \frac{dw}{dt} = \frac{Q'\,e^{i\theta}}{2w^2},
```

so the branch of the square root is transported by the ODE and never chosen discretely. The
combinatorial content is read off with [`topology_signature`](@ref), [`finite_lines`](@ref),
[`edges`](@ref) and [`n_infinite_lines`](@ref), and that signature is what the
[cluster bridge](bridge.md) consumes.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/stokes_graph.jl"]
```

## Saddles

At special angles a Stokes line runs from one turning point into another: a *saddle
connection*, and the topology jumps there. These are the BPS states: the critical angle is the
phase of the central charge ``Z_\gamma = 2\int\sqrt{Q}``, and ``|Z_\gamma|`` is its mass.
[`saddles`](@ref) locates them by bisecting in ``\theta`` between candidate angles, and
[`stokes_graph_family`](@ref) sweeps a range of angles so you can watch the topology flip
across one.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/saddles.jl"]
```
