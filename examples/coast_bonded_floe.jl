"""
One square floe bonded to a wall (a topography element) along their shared face, pushed
away from the wall at 5 m/s. Coupling is off, so the bond is the only force pulling it back.
Expected Behavior:
    Floe moves away from the wall, the bond stretches and pulls it back, and the floe
    oscillates along x (bouncing off the wall when it reaches it) while bond damping
    slowly removes energy.
"""

using Subzero, StructArrays, StaticArrays
using CairoMakie, GeoInterfaceMakie   # loads Subzero's plotting extension for plot_sim

const FT = Float64
const side = 1e4      # floe side length (m)
const height = 1.0    # floe thickness (m)
const Δt = 10
const nΔt = 2000

# domain sized so the floe and wall fill a reasonable part of the video frame
grid = RegRectilinearGrid(; x0 = -1e4, xf = 3e4, y0 = -1.5e4, yf = 2.5e4, Δx = 1e4, Δy = 1e4)
# wall is a 2 km thick strip whose right face is at x = 0
wall = [[[[-2e3, -1e4], [-2e3, 2e4], [0.0, 2e4], [0.0, -1e4], [-2e3, -1e4]]]]
domain = Subzero.Domain(;
    north = OpenBoundary(North; grid),
    south = OpenBoundary(South; grid),
    east = OpenBoundary(East; grid),
    west = OpenBoundary(West; grid),
    topography = initialize_topography_field(; coords = wall),
)
ocean = Ocean(; grid, u = 0.0, v = 0.0, temp = 0.0)
atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = 0.0)

# floe sits right of the wall, touching it along x = 0, moving away from it
floe_settings = FloeSettings()
square = [[[0.0, 0.0], [0.0, side], [side, side], [side, 0.0], [0.0, 0.0]]]
F = Floe{FT}(square, height; floe_settings, id = 1, u = 5.0)

# bond anchored at the middle of the shared face, normal pointing from the floe to the wall
Subzero.add_coast_bond!(F, SVector(0.0, side / 2), SVector(-1.0, 0.0), side)

fracture_settings = FractureSettings(;
    fractures_on = true,
    bond_criteria = MohrCoulombBondCriteria(;
        floe  = MohrCoulomb(; cohesion = 5.0e4, μ = 0.92),
        coast = MohrCoulomb(; cohesion = 2e4, μ = 0.5, tensile_strength = 1e5),
    ),
)

model = Model(; grid, ocean, atmos, domain, floes = StructArray([F]))

dir = joinpath(@__DIR__, "output", "coast_bonded_floe")
writers = OutputWriters(
    InitialStateOutputWriter(; dir, overwrite = true),
    FloeOutputWriter(10; dir, overwrite = true),
)
simulation = Simulation(;
    name = "coast_bonded_floe",
    model,
    Δt,
    nΔt,
    floe_settings,
    coupling_settings = CouplingSettings(coupling_on = false),
    bond_settings = BondSettings(bonds_on = true, coastline_bonds = true),
    fracture_settings,
    writers,
)

run!(simulation)

plot_sim(
    joinpath(dir, "floes.jld2"),
    joinpath(dir, "initial_state.jld2"),
    Δt,
    joinpath(dir, "coast_bonded_floe.mp4"),
)
