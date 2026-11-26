"""
Two square floes bonded along their shared face, pushed apart at 5 m/s each.
Coupling is off, so the bond is the only force pulling them back together.
Expected Behavior:
    Floes separate, the bond stretches and pulls them back, and the motion
    oscillates along x while bond damping slowly removes energy.
"""

using Subzero, StructArrays, StaticArrays
using CairoMakie, GeoInterfaceMakie   # loads Subzero's plotting extension for plot_sim

const FT = Float64
const side = 1e4      # floe side length (m)
const height = 1.0    # floe thickness (m)
const Δt = 10
const nΔt = 2000

# domain sized so the floes fill a reasonable part of the video frame
grid = RegRectilinearGrid(; x0 = -3e4, xf = 3e4, y0 = -2.5e4, yf = 3.5e4, Δx = 1e4, Δy = 1e4)
domain = Subzero.Domain(;
    north = OpenBoundary(North; grid),
    south = OpenBoundary(South; grid),
    east = OpenBoundary(East; grid),
    west = OpenBoundary(West; grid),
)
ocean = Ocean(; grid, u = 0.0, v = 0.0, temp = 0.0)
atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = 0.0)

# A sits left of B; they share the face x = 0
square(x0) = [[[x0, 0.0], [x0, side], [x0 + side, side], [x0 + side, 0.0], [x0, 0.0]]]
floe_settings = FloeSettings()
A = Floe{FT}(square(-side), height; floe_settings, id = 1, u = -5.0)
B = Floe{FT}(square(0.0), height; floe_settings, id = 2, u = 5.0)

# bond anchored at the middle of the shared face, normal pointing from A to B
Subzero.add_bond!(A, B, SVector(0.0, side / 2), SVector(1.0, 0.0), side)

fracture_settings = FractureSettings(;
    fractures_on = true,
    bond_criteria = MohrCoulombBondCriteria(;
        floe  = MohrCoulomb(; cohesion = 5.5e4, μ = 0.92),
        coast = MohrCoulomb(; cohesion = 2e5, μ = 0.5, tensile_strength = 1e5),
    ),
)

model = Model(; grid, ocean, atmos, domain, floes = StructArray([A, B]))

dir = joinpath(@__DIR__, "output", "bonded_floes")
writers = OutputWriters(
    InitialStateOutputWriter(; dir, overwrite = true),
    FloeOutputWriter(10; dir, overwrite = true),
)
simulation = Simulation(;
    name = "bonded_floes",
    model,
    Δt,
    nΔt,
    floe_settings,
    coupling_settings = CouplingSettings(coupling_on = false),
    bond_settings = BondSettings(bonds_on = true),
    fracture_settings,
    writers,
)

run!(simulation)

plot_sim(
    joinpath(dir, "floes.jld2"),
    joinpath(dir, "initial_state.jld2"),
    Δt,
    joinpath(dir, "bonded_floes.mp4"),
)
