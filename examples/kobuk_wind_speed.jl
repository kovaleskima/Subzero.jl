"""
Wind speed sweep of the Kobuk Lake simulation (see kobuk_lake_setup.jl). The simulation is
run 8 times with wind blowing toward the east at 2, 4, ..., 16 m/s, on the same floe field.
For each run, the bond failure rate is the largest number of bonds broken at any point in the
run divided by the number of bonds when the run started (floe-floe and coast bonds together;
each floe-floe bond counts once). Rates are saved to failure_rates.csv as each run finishes,
the number of broken bonds through each run to broken_<speed>ms.csv, and a plot of failure
rate against wind speed to failure_rate_vs_speed.png, all in output/kobuk_wind_speed.
Each run takes as long as kobuk_lake.jl; set the environment variable KOBUK_NSTEPS to run
fewer timesteps, e.g. `KOBUK_NSTEPS=200 julia --project=. -t auto examples/kobuk_wind_speed.jl`.
Expected Behavior:
    Stronger winds push the ice cover harder against the eastern shore, so the failure rate
    should increase with wind speed.
"""

using CairoMakie
include("kobuk_lake_setup.jl")

const wind_speeds = 2.0:2.0:16.0    # m/s, wind blowing toward the east
const nsteps = parse(Int, get(ENV, "KOBUK_NSTEPS", "21600"))

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

dir = joinpath(@__DIR__, "output", "kobuk_wind_speed")
mkpath(dir)
rates_file = joinpath(dir, "failure_rates.csv")
open(rates_file, "w") do io
    println(io, "wind_speed,u_wind,v_wind,initial_floe_bonds,initial_coast_bonds,",
        "max_broken_floe_bonds,max_broken_coast_bonds,max_broken_bonds,failure_rate")
end

rates = Float64[]
for speed in wind_speeds
    u_wind, v_wind = speed, 0.0
    label = "$(round(Int, speed))ms"
    println("Wind toward E at $speed m/s")
    sim = kobuk_lake_simulation(; u_wind, v_wind, nΔt = nsteps,
        name = "kobuk_wind_$label", verbose = false)
    nfloe0, ncoast0, broken = run_counting_bonds!(sim)
    total = vec(sum(broken; dims = 2))
    max_broken = maximum(total)
    rate = max_broken / (nfloe0 + ncoast0)
    push!(rates, rate)
    println("  $max_broken of $(nfloe0 + ncoast0) bonds broken, failure rate = $(round(rate; digits = 4))")
    open(rates_file, "a") do io
        println(io, join((speed, u_wind, v_wind, nfloe0, ncoast0, maximum(broken[:, 1]),
            maximum(broken[:, 2]), max_broken, rate), ","))
    end
    open(joinpath(dir, "broken_$label.csv"), "w") do io
        println(io, "tstep,broken_floe_bonds,broken_coast_bonds")
        for k in axes(broken, 1)
            println(io, k - 1, ",", broken[k, 1], ",", broken[k, 2])
        end
    end
end

# failure rate against wind speed
pct = 100 .* rates
fig = Figure(size = (700, 500))
ax = Axis(fig[1, 1];
    title = "Bond failure rate by wind speed\nwind toward east, $(round(nsteps * Δt / 3600; digits = 2)) h simulated",
    xlabel = "Wind speed (m/s)",
    ylabel = "Bond failure (%)",
    xticks = collect(wind_speeds),
)
ylims!(ax, 0, nothing)
lines!(ax, collect(wind_speeds), pct; color = :steelblue, linewidth = 2)
scatter!(ax, collect(wind_speeds), pct; color = :steelblue, markersize = 12)
save(joinpath(dir, "failure_rate_vs_speed.png"), fig)
println("Saved results to $dir")
