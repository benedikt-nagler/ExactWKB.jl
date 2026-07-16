# Makie weakdep extension: Stokes-graph plotting.
#
# Renders a StokesGraph in the z-plane - turning points as markers, traced Stokes lines
# as curves (saddles highlighted), with a DataAspect so angles read true. A vector
# method lays out a θ-family as a grid of panels.
module ExactWKBMakieExt

using ExactWKB
using ExactWKB: lines, points, location, is_finite_line, turning_points, is_simple
using Makie: Makie, Figure, Axis, lines!, scatter!, text!, DataAspect, axislegend, @L_str

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

const _MARKED_COLOR = :darkorange

function ExactWKB.plot_triangulation(g::StokesGraph, t::IdealTriangulation;
                                     resolution = (600, 600), title = nothing,
                                     legend = true)
    fig = Figure(size = resolution)
    ttl = title === nothing ?
          "Stokes graph + dual triangulation  (θ = $(round(Float64(g.theta), digits = 4)))" :
          title
    ax = Axis(fig[1, 1]; title = ttl, xlabel = L"\mathrm{Re}\,z",
              ylabel = L"\mathrm{Im}\,z", aspect = DataAspect())
    _draw!(ax, g)
    # the boundary circle (where the rays exit) with the marked points at the
    # asymptotic directions, and the diagonals as chords
    R = maximum(abs(points(l)[end]) for l in ExactWKB.infinite_lines(g))
    φ = range(0, 2π; length = 256)
    lines!(ax, R .* cos.(φ), R .* sin.(φ); color = (:gray, 0.6), linestyle = :dot)
    mp = [R * cis(a) for a in t.marked_angles]
    scatter!(ax, Float64.(real.(mp)), Float64.(imag.(mp));
             color = _MARKED_COLOR, markersize = 12, label = "marked point")
    first_diag = true
    for e in eachindex(t.edge_endpoints)
        t.is_diagonal[e] || continue
        zu, zv = mp[t.edge_endpoints[e][1]], mp[t.edge_endpoints[e][2]]
        lines!(ax, Float64[real(zu), real(zv)], Float64[imag(zu), imag(zv)];
               color = _MARKED_COLOR, linewidth = 2, linestyle = :dash,
               label = first_diag ? "diagonal" : nothing)
        first_diag = false
        mid = (zu + zv) / 2
        (i, j) = t.diagonal_tp_pair[e]
        text!(ax, Float64(real(mid)), Float64(imag(mid));
              text = "γ($i,$j)", color = _MARKED_COLOR, fontsize = 12)
    end
    legend && axislegend(ax; position = :rt, framevisible = false)
    fig
end

end # module ExactWKBMakieExt
