"""
Uniaxial compression of a bonded pack of ~1000 Voronoi-tessellated floes filling a 20 km by
20 km box. The northern boundary moves south, the southern boundary is a fixed wall, and the
domain extends 5 km past the pack on the east and west, so the pack is compressed along y and
free to bow outward along x before reaching the eastern/western walls. With bonds on, every pair of floes sharing a face is bonded along it
when the simulation starts (no coast bonds). Bonds break with Mohr-Coulomb bond criteria;
floes themselves don't fracture. The ocean and atmosphere are at rest, so coupling only adds
drag.
Expected Behavior:
    The bonded pack first compresses elastically and bows outward into the open water on its
    sides, then bonds begin to fail and fractures propagate through the pack, breaking it into
    bonded aggregates that jam against each other and slide out sideways.
"""

using Subzero, Random, Statistics
using CairoMakie, GeoInterfaceMakie   # loads Subzero's plotting extension for plot_sim

const FT = Float64
const Lx = 2e4          # pack width (m)
const Ly = 2e4          # pack and domain height (m)
const pad = 5e3         # open water between the pack and the eastern/western walls (m)
const Δgrid = 1e3       # ocean/atmos grid spacing (m)
const nfloes = 1000
const hmean = 1.0       # mean floe thickness (m)
const Δh = 0.25         # floe thickness range around the mean (m)
const v_wall = 0.05      # speed of the northern boundary (m/s)
const Δt = 2            # larger timesteps make small, stiffly bonded floes rattle
const nΔt = 10000       # northern boundary moves 2 km, i.e. 10% strain

grid = RegRectilinearGrid(; x0 = -pad, xf = Lx + pad, y0 = 0.0, yf = Ly, Δx = Δgrid, Δy = Δgrid)
domain = Subzero.Domain(;
    north = MovingBoundary(North; grid, v = -v_wall),
    south = CollisionBoundary(South; grid),
    east = CollisionBoundary(East; grid),
    west = CollisionBoundary(West; grid),
)
ocean = Ocean(; grid, u = 0.0, v = 0.0, temp = 0.0)
atmos = Atmos(; grid, u = 0.0, v = 0.0, temp = 0.0)

# floes average ~0.4 km², so lower the minimum floe area from its 1 km² default
floe_settings = FloeSettings(; min_floe_area = 1e4)
generator = VoronoiTesselationFieldGenerator(; nfloes, concentrations = [1.0], hmean, Δh)
# tessellate only the pack, leaving the padding on either side as open water
floe_bounds = Subzero.make_polygon([[[0.0, 0.0], [0.0, Ly], [Lx, Ly], [Lx, 0.0], [0.0, 0.0]]], FT)
floes = initialize_floe_field(; generator, domain, floe_bounds, floe_settings, rng = Xoshiro(1))

# Young's modulus scaled to floe size (see Constants)
consts = Constants(; E = 2.5e3 * (mean(sqrt.(floes.area)) + minimum(sqrt.(floes.area))))

fracture_settings = FractureSettings(;
    fractures_on = true,
    bond_criteria = MohrCoulombBondCriteria(; floe = MohrCoulomb(; cohesion = 2.0e4, μ = 0.92)),
)

model = Model(; grid, ocean, atmos, domain, floes)

dir = joinpath(@__DIR__, "output", "uniaxial_compression")
writers = OutputWriters(
    InitialStateOutputWriter(; dir, overwrite = true),
    FloeOutputWriter(100; dir, overwrite = true),
)
simulation = Simulation(;
    name = "uniaxial_compression",
    model,
    consts,
    Δt,
    nΔt,
    verbose = true,
    rng = Xoshiro(1),
    floe_settings,
    # bond stiffness = contact stiffness / ratio, separately across the face (tension and
    # compression) and along it (shear)
    bond_settings = BondSettings(;
        bonds_on = true,
        normal_stiffness_ratio = 1,
        shear_stiffness_ratio = 1,
    ),
    fracture_settings,
    writers,
)

run!(simulation)

plot_sim(
    joinpath(dir, "floes.jld2"),
    joinpath(dir, "initial_state.jld2"),
    Δt,
    joinpath(dir, "uniaxial_compression.mp4"),
)
