"""
Wind-driven breakup of a bonded ice cover on Kobuk Lake. The lake's coastline is topography
and the ice initially fills the lake region from Kobuk.mat (see data/kobuk/export_kobuk.m),
tessellated into ~1000 Voronoi floes. Every pair of floes sharing a face is bonded along it
when the simulation starts, and floes along the shore are bonded to the coast. Bonds break
with Mohr-Coulomb bond criteria; floes themselves don't fracture. A constant, uniform 10 m/s
wind blows west across the domain over an ocean at rest, which is the only forcing. The
setup is in kobuk_lake_setup.jl. The animation draws the floes and bonds over a satellite image
of the lake (data/kobuk/kobuk_satellite.jpg, made by data/kobuk/make_kobuk_satellite.py).
Expected Behavior:
    The wind pushes the ice cover west, but the lake is almost fully ice covered, so after a
    brief drift the cover jams against the western shore and barely moves. With these bond
    strengths only a few bonds break, at narrow necks of ice near the western shore; the
    rest of the cover stays intact for the 6 hours simulated.
"""

using CairoMakie, GeoInterfaceMakie   # loads Subzero's plotting extension for plot_sim
include("kobuk_lake_setup.jl")

dir = joinpath(@__DIR__, "output", "kobuk_lake")
writers = OutputWriters(
    InitialStateOutputWriter(; dir, overwrite = true),
    FloeOutputWriter(100; dir, overwrite = true),
)
simulation = kobuk_lake_simulation(; u_wind = 10, v_wind = 0.0, nΔt = 43200, writers)

run!(simulation)

plot_sim(
    joinpath(dir, "floes.jld2"),
    joinpath(dir, "initial_state.jld2"),
    Δt,
    joinpath(dir, "kobuk_lake.mp4");
    background = Makie.FileIO.load(joinpath(kobuk_data_dir, "kobuk_satellite.jpg")),
)
