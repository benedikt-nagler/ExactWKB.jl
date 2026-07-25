```@meta
CurrentModule = ExactWKB
```

# Errors

Every invalid operation throws a typed error from the [`ExactWKBError`](@ref) hierarchy, each
with an informative `showerror`. Generic bad keyword arguments reuse
`Resurgence.InvalidArgument`.

The design rule behind them: a computation that has left the regime where it is valid raises
rather than returns. A graph with a saddle connection is not a generic graph, a coalescing pair
of turning points is not a spectral cycle, and a non-decaying wall symbol is not a wall symbol -
in each case the alternative would be a plausible-looking number that means nothing.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/errors.jl"]
```
