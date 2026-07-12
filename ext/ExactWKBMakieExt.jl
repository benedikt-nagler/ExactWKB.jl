# Makie weakdep extension: Stokes-graph plotting.
#
# Renders a StokesGraph in the z-plane — turning points as markers, traced Stokes lines
# as curves (saddles highlighted), with a DataAspect so angles read true. A vector
# method lays out a θ-family as a grid of panels.
module ExactWKBMakieExt

using ExactWKB
using ExactWKB: lines, points, location, is_finite_line, turning_points, is_simple
using Makie: Makie, Figure, Axis, lines!, scatter!, DataAspect, axislegend, @L_str

# colour roles (backend-agnostic named colours)
const _INF_COLOR = :steelblue
const _SADDLE_COLOR = :crimson
const _TP_COLOR = :black

function _draw!(ax, g::StokesGraph)
    saddle_drawn = false
    inf_drawn = false
    for l in lines(g)
        pts = points(l)
        xs = Float64.(real.(pts)); ys = Float64.(imag.(pts))
        if is_finite_line(l)
            lines!(ax, xs, ys; color = _SADDLE_COLOR, linewidth = 2.5,
                   label = saddle_drawn ? nothing : "saddle")
            saddle_drawn = true
        else
            lines!(ax, xs, ys; color = _INF_COLOR, linewidth = 1.2,
                   label = inf_drawn ? nothing : "→ ∞")
            inf_drawn = true
        end
    end
    tp = turning_points(g)
    scatter!(ax, [Float64(real(location(t))) for t in tp],
             [Float64(imag(location(t))) for t in tp];
             color = _TP_COLOR, marker = :xcross, markersize = 14, label = "turning point")
    ax
end

function ExactWKB.plot_stokes_graph(g::StokesGraph; resolution = (600, 600),
                                    title = nothing, legend = true)
    fig = Figure(size = resolution)
    ttl = title === nothing ? "Stokes graph  (θ = $(round(Float64(g.theta), digits = 4)))" : title
    ax = Axis(fig[1, 1]; title = ttl, xlabel = L"\mathrm{Re}\,z", ylabel = L"\mathrm{Im}\,z",
              aspect = DataAspect())
    _draw!(ax, g)
    legend && axislegend(ax; position = :rt, framevisible = false)
    fig
end

function ExactWKB.plot_stokes_graph(gs::AbstractVector{<:StokesGraph};
                                    resolution = nothing, ncols = nothing)
    n = length(gs)
    cols = ncols === nothing ? ceil(Int, sqrt(n)) : ncols
    rows = ceil(Int, n / cols)
    res = resolution === nothing ? (320 * cols, 320 * rows) : resolution
    fig = Figure(size = res)
    for (k, g) in enumerate(gs)
        r, c = fldmod1(k, cols)
        ax = Axis(fig[r, c]; title = "θ = $(round(Float64(g.theta), digits = 4))",
                  aspect = DataAspect())
        _draw!(ax, g)
    end
    fig
end

end # module ExactWKBMakieExt
