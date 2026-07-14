# ExactWKB.jl

Exact WKB analysis for one-dimensional Schrödinger equations `ħ²ψ″ = Q(z)ψ` with polynomial `Q`. The package computes turning points, the all-orders WKB expansion, quantum periods and Voros symbols, and Stokes graphs with their saddle connections.

Borel summation of the resulting divergent series is provided by
[Resurgence.jl](../Resurgence.jl).

```julia
using ExactWKB

prob = SchrodingerProblem([3//4, 0, -2, 0, 1])   # double well Q = (z²−1)² − 1/4
tps  = turning_points(prob)
w    = wkb_expansion(prob; order = 6)
vs   = voros_symbol(w, encircling_contour(tps[2], tps[3]))
g    = stokes_graph(prob; theta = 0.0)
```

`ExactWKB.jl` depends on the unregistered packages `Resurgence.jl` and
`ClusterAlgebras.jl`; install all three with `Pkg.develop` from a local checkout.
See `examples/m3_flagship.jl` for a complete example. MIT license.
