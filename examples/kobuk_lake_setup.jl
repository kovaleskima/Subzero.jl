#= Kobuk Lake setup shared by kobuk_lake.jl and kobuk_wind_direction.jl: the lake's coastline
is topography and the ice initially fills the lake region from Kobuk.mat (see
data/kobuk/export_kobuk.m), tessellated into ~1000 Voronoi floes. Every pair of floes sharing
a face is bonded along it when the simulation starts, and floes along the shore are bonded to
the coast. Bonds break with Mohr-Coulomb bond criteria; floes themselves don't fracture. A
constant, uniform wind over an ocean at rest is the only forcing. =#

using Subzero, Random, Statistics

const FT = Float64
const Δgrid = 1e3       # ocean/atmos grid spacing (m)
const nfloes = 1000
const hmean = 0.25      # floe thickness (m)
const Δh = 0.0          # floe thickness range around the mean (m)
const Δt = 2            # larger timesteps make small, stiffly bonded floes rattle

# Polygons from a CSV with rows `polygon id, x, y`, as coordinates of closed rings
function read_polygons(fn)
    rings = Dict{Int, Vector{Vector{FT}}}()
    for line in eachline(fn)
        id, x, y = parse.(FT, split(line, ','))
        push!(get!(rings, Int(id), Vector{FT}[]), [x, y])
    end
    return [[push!(rings[id], rings[id][1])] for id in sort!(collect(keys(rings)))]
end
const kobuk_data_dir = joinpath(@__DIR__, "data", "kobuk")
const coast_coords = read_polygons(joinpath(kobuk_data_dir, "kobuk_coast.csv"))  # land around the lake
const ice_coords = read_polygons(joinpath(kobuk_data_dir, "kobuk_ice.csv"))      # initial ice region

"""
    kobuk_lake_simulation(; u_wind, v_wind, nΔt, name, writers, verbose)

Kobuk Lake simulation with a constant, uniform wind of velocity (`u_wind`, `v_wind`) (m/s).
The floe field is the same for every call.
"""
function kobuk_lake_simulation(;
    u_wind = 0.0,
    v_wind = 10.0,
    nΔt = 43200,
    name = "kobuk_lake",
    writers = OutputWriters(),
    verbose = true,
)
    # the coastline polygons fill an 80 km by 70 km box around the lake
    grid = RegRectilinearGrid(; x0 = -3.5e4, xf = 3.5e4, y0 = -3.5e4, yf = 3.5e4, Δx = Δgrid, Δy = Δgrid)
    domain = Subzero.Domain(;
        north = CollisionBoundary(North; grid),
        south = CollisionBoundary(South; grid),
        east = CollisionBoundary(East; grid),
        west = CollisionBoundary(West; grid),
        topography = initialize_topography_field(FT; coords = coast_coords),
    )
    ocean = Ocean(; grid, u = 0.0, v = 0.0, temp = 0.0)
    atmos = Atmos(; grid, u = u_wind, v = v_wind, temp = 0.0)

    # floes average ~1 km², so lower the minimum floe area from its 1 km² default
    floe_settings = FloeSettings(; min_floe_area = 1e3)
    generator = VoronoiTesselationFieldGenerator(; nfloes, concentrations = [1.0], hmean, Δh)
    # tessellate only the initial ice region
    floe_bounds = Subzero.make_multipolygon(ice_coords, FT)
    floes = initialize_floe_field(; generator, domain, floe_bounds, floe_settings, rng = Xoshiro(1))
    # the ice region and coastline were outlined separately, so shore floes sit up to a few
    # hundred meters from the coast; bond faces within 400 m of it, which leaves out faces
    # bordering the open water between the ice and the eastern shore
    #ncoast = Subzero.bond_floes_to_coast!(floes, domain; max_gap = 400.0)
    #verbose && println("Bonded $ncoast floe faces to the coast.")

    # Young's modulus scaled to floe size (see Constants)
    consts = Constants(; E = 2.5e3 * (mean(sqrt.(floes.area)) + minimum(sqrt.(floes.area))))

    fracture_settings = FractureSettings(;
        fractures_on = true,
        bond_criteria = MohrCoulombBondCriteria(; floe = MohrCoulomb(; cohesion = 0.8e4, μ = 0.92)),
    )

    model = Model(; grid, ocean, atmos, domain, floes)
    return Simulation(;
        name,
        model,
        consts,
        Δt,
        nΔt,
        verbose,
        rng = Xoshiro(1),
        floe_settings,
        # bond stiffness = contact stiffness / ratio, separately across the face (tension and
        # compression) and along it (shear)
        bond_settings = BondSettings(;
            bonds_on = true,
            coastline_bonds = false,
            normal_stiffness_ratio = 2,
            shear_stiffness_ratio = 1,
        ),
        fracture_settings,
        writers,
    )
end
