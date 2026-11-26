"""
Wind direction sweep of the Kobuk Lake simulation (see kobuk_lake_setup.jl). The simulation
is run 8 times with a 10 m/s wind toward each cardinal and ordinal direction, on the same
floe field. For each run, the bond failure rate is the largest number of bonds broken at any
point in the run divided by the number of bonds when the run started (floe-floe and coast
bonds together; each floe-floe bond counts once). Rates are saved to failure_rates.csv as
each run finishes, the number of broken bonds through each run to broken_<direction>.csv, and
a polar plot of failure rate against wind direction to failure_rate_polar.png, all in
output/kobuk_wind_direction.
Each run takes as long as kobuk_lake.jl; set the environment variable KOBUK_NSTEPS to run
fewer timesteps, e.g. `KOBUK_NSTEPS=200 julia --project=. -t auto examples/kobuk_wind_direction.jl`.
Expected Behavior:
    Each run pushes the ice cover toward a different shore, so the failure rate varies with
    wind direction. Positive u and v push the ice east and north.
"""

using CairoMakie
include("kobuk_lake_setup.jl")

const wind_speed = 10.0     # m/s
const nsteps = parse(Int, get(ENV, "KOBUK_NSTEPS", "10800"))
# directions the wind blows toward, as angles counterclockwise from east (so north is 90°)
const directions = ["E" => 0, "NE" => 45, "N" => 90, "NW" => 135,
                    "W" => 180, "SW" => 225, "S" => 270, "SE" => 315]

# number of floe-floe bonds (each bond is stored on both of its floes) and coast bonds
count_bonds(floes) = (sum(b -> length(b.floe), floes.bonds) ÷ 2, sum(b -> length(b.coast), floes.bonds))

# Run the simulation like run! does, counting bonds after every timestep. Returns the number
# of floe-floe and coast bonds at the start and the number broken through each timestep.
function run_counting_bonds!(sim)
    Subzero.startup_sim(sim)    # bonds adjacent floes
    nfloe0, ncoast0 = count_bonds(sim.model.floes)
    broken = zeros(Int, sim.nΔt + 1, 2)
    for tstep in 0:sim.nΔt
        Subzero.timestep_sim!(sim, tstep)
        nfloe, ncoast = count_bonds(sim.model.floes)
        broken[tstep + 1, :] .= (nfloe0 - nfloe, ncoast0 - ncoast)
    end
    Subzero.teardown_sim(sim)
    return nfloe0, ncoast0, broken
end

dir = joinpath(@__DIR__, "output", "kobuk_wind_direction")
mkpath(dir)
rates_file = joinpath(dir, "failure_rates.csv")
open(rates_file, "w") do io
    println(io, "direction,angle_deg,u_wind,v_wind,initial_floe_bonds,initial_coast_bonds,",
        "max_broken_floe_bonds,max_broken_coast_bonds,max_broken_bonds,failure_rate")
end

rates = Float64[]
for (name, θ) in directions
    u_wind, v_wind = wind_speed .* (cosd(θ), sind(θ))
    println("Wind toward $name ($(θ)°): u = $(round(u_wind; digits = 2)), v = $(round(v_wind; digits = 2)) m/s")
    sim = kobuk_lake_simulation(; u_wind, v_wind, nΔt = nsteps,
        name = "kobuk_wind_$name", verbose = false)
    nfloe0, ncoast0, broken = run_counting_bonds!(sim)
    total = vec(sum(broken; dims = 2))
    max_broken = maximum(total)
    rate = max_broken / (nfloe0 + ncoast0)
    push!(rates, rate)
    println("  $max_broken of $(nfloe0 + ncoast0) bonds broken, failure rate = $(round(rate; digits = 4))")
    open(rates_file, "a") do io
        println(io, join((name, θ, u_wind, v_wind, nfloe0, ncoast0, maximum(broken[:, 1]),
            maximum(broken[:, 2]), max_broken, rate), ","))
    end
    open(joinpath(dir, "broken_$name.csv"), "w") do io
        println(io, "tstep,broken_floe_bonds,broken_coast_bonds")
        for k in axes(broken, 1)
            println(io, k - 1, ",", broken[k, 1], ",", broken[k, 2])
        end
    end
end

# polar plot of failure rate against wind direction; angles are counterclockwise from east,
# so north is at the top as on a map
θs = [deg2rad(θ) for (_, θ) in directions]
pct = 100 .* rates
fig = Figure(size = (650, 650))
ax = PolarAxis(fig[1, 1];
    title = "Bond failure rate (%) by wind direction\n$(wind_speed) m/s wind toward each direction, $(round(nsteps * Δt / 3600; digits = 2)) h simulated",
    thetaticks = (θs, first.(directions)),
    rlimits = (0, nothing),
    rticks = Makie.WilkinsonTicks(4),
    rtickangle = deg2rad(22.5),   # radial tick labels between the E and NE spokes
)
lines!(ax, [θs; θs[1]], [pct; pct[1]]; color = :steelblue, linewidth = 2)
scatter!(ax, θs, pct; color = :steelblue, markersize = 12)
save(joinpath(dir, "failure_rate_polar.png"), fig)
println("Saved results to $dir")
