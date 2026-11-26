module SubzeroMakieExt

using CairoMakie, GeoInterfaceMakie
using Subzero, JLD2
import Subzero: prettytime, plot_sim

"""
    plot_sim(floe_fn, initial_state_fn, Δt, output_fn; background = nothing)

Basic plotting of a simulation using the simulation's floe and initial state files. This
function is meant for basic plotting during testing and as an example of how to create
a video using Makie. The user can (and should!) write their own plotting code to their own
specifications and needs.

This code does not have an underlying ocean, but there are comments on how to add an ocean
heatmap and colorbar.

## Arguments:
- `floe_fn::String`: $(Subzero.FLOE_FN_DEF)
- `initial_state_fn::String`: $(Subzero.INITIAL_STATE_FN_DEF)
- $(Subzero.ΔT_DEF)
- `output_fn::String`: $(Subzero.MP4_OUTPUT_FN)

## Keyword arguments:
- `background`: an image (matrix of colors, first row at the top) covering the domain from
    its west to east and south to north walls, drawn under the floes in place of the
    topography, e.g. a satellite image of the region.
"""
function plot_sim(
    floe_fn,
    initial_state_fn,
    Δt,
    output_fn;
    background = nothing,
)
    # Color Information
    topo_color = RGBf(147/255, 129/255, 102/255)
    floe_color = RGBf(217/255, 226/255, 225/255)
    ocean_color = RGBf(4/255, 31/255, 74/255)
    # Domain Information
    domain = load(initial_state_fn)["sim"].model.domain
    xmax, xmin = domain.east.val, domain.west.val
    ymax, ymin = domain.north.val, domain.south.val
    Δx, Δy = xmax - xmin, ymax - ymin
    # Floe Information
    file = jldopen(floe_fn)
    sim_polys = file["poly"]
    timesteps = keys(sim_polys)
    # Set up observables for recording (updated whenever `time[]` is set to a new value)
    time = Observable(timesteps[1])  # note these are strings as we index into a JLD2 file
    time_polys = @lift(sim_polys[$time])
    title_string = @lift(Subzero.prettytime(parse(Float64, $time) * Δt))
    # Bond Information (only if bonds were output)
    plot_bonds = all(k -> haskey(file, k), ("bonds", "centroid", "α", "id", "ghost_id"))
    if plot_bonds
        unbroken_pts, broken_pts = bond_points(file, timesteps)
        time_unbroken = @lift(unbroken_pts[$time])
        time_broken = @lift(broken_pts[$time])
    end

    #=
    Note that if you had Nx by Ny by time ocean data, this would be a good place to connect
    it to the time observable! To learn more about Observables and how to use them with
    animations see here: https://docs.makie.org/stable/explanations/animation
    =#

    # Set up figure
    fig = Figure()
    Axis(
        fig[1, 1];
        backgroundcolor = ocean_color,
        limits = (xmin, xmax, ymin, ymax),
        title = title_string,
        xlabel = "Meters",
        xticklabelrotation = pi/4,
        ylabel = "Meters",
        yticklabelrotation = pi/4,
    )

    # Plot starting state (background + floes + topography)
    if !isnothing(background)  # image rows run top to bottom, Makie's x index runs west to east
        image!((xmin, xmax), (ymin, ymax), rotr90(background))
    end
    poly!(time_polys, color = floe_color, strokecolor = floe_color, strokewidth = 0.5)  # floes
    if !isempty(domain.topography) && isnothing(background)  # topography
        topo_mp = domain.topography.poly
        poly!(topo_mp, color = topo_color)
    end
    if plot_bonds  # unbroken bonds as red dots, broken bonds as blue x's
        scatter!(time_unbroken, color = :red, marker = :circle, markersize = 2)
        scatter!(time_broken, color = :blue, marker = :xcross, markersize = 3)
    end

    #=
    Note that if you had Nx by Ny by time ocean data, this would be a good place to plot it
    with a heat map and to add a color bar on the fig[1,2] axis! Check out Makie
    documentation here: https://docs.makie.org/stable/reference/plots/heatmap 
    =#

    # Resize Figure
    colsize!(fig.layout, 1, Aspect(1, Δx / Δy))
    resize_to_layout!(fig)

    # Create movie
    record(fig, output_fn, timesteps; framerate = 20) do t
        time[] = t  # updating the time auto updates the related time_polys
    end

    # End movie and close file
    close(file)
end

# world-frame location of a bond anchor given the floe's centroid, rotation, and body-frame
# attachment vector
_bond_world(c, α, a) = Point2f(c[1] + cos(α) * a[1] - sin(α) * a[2],
    c[2] + sin(α) * a[1] + cos(α) * a[2])

#=
For each output timestep, find the locations of unbroken and broken bonds. A bond is broken
at a timestep if it existed at an earlier timestep but no longer does. Broken bonds keep
moving with the floe they were attached to; if that floe no longer exists, they stay at
their last location.
=#
function bond_points(file, timesteps)
    unbroken = Dict{String, Vector{Point2f}}()
    broken = Dict{String, Vector{Point2f}}()
    # bond key => (floe id, body-frame attachment, last world location)
    seen = Dict{Any, Tuple{Int, Any, Point2f}}()
    for t in timesteps
        ids, ghost_ids = file["id"][t], file["ghost_id"][t]
        centroids, αs, bonds = file["centroid"][t], file["α"][t], file["bonds"][t]
        idx = Dict(ids[i] => i for i in eachindex(ids) if ghost_ids[i] == 0)
        alive = Set{Any}()
        unbroken[t] = Point2f[]
        for (id, i) in idx
            for b in bonds[i].floe
                id < b.partner_id || continue  # each floe-floe bond is stored on both floes
                key = (id, b.partner_id)
                p = _bond_world(centroids[i], αs[i], b.attach_self)
                push!(alive, key); push!(unbroken[t], p)
                seen[key] = (id, b.attach_self, p)
            end
            for b in bonds[i].coast
                key = (id, Tuple(b.anchor))
                p = Point2f(b.anchor[1], b.anchor[2])
                push!(alive, key); push!(unbroken[t], p)
                seen[key] = (id, nothing, p)
            end
        end
        broken[t] = Point2f[]
        for (key, (id, a, p)) in seen
            key in alive && continue
            i = get(idx, id, 0)
            push!(broken[t], (i == 0 || isnothing(a)) ? p : _bond_world(centroids[i], αs[i], a))
        end
    end
    return unbroken, broken
end

end  # module
