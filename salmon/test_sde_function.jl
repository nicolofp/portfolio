using Statistics, Distributions, LinearAlgebra, Random, Plots

function simulate_popfish(N0::Int, r::Float64, K::Float64, sigma::Float64, 
                          Wh::Float64, mu::Float64, T::Float64, dt::Float64)
    rng = MersenneTwister(90)
    times = collect(0.0:dt:T)
    X = rand(rng, Uniform(0.5, 1.0), N0)
    remaining = trues(N0)
    results = [(0.0, N0, mean(X), std(X), 0, 0.0, 0.0)]
    harvest_times = Float64[]
    death_times = Float64[]

    for t_idx in 2:length(times)
        t = times[t_idx]
        if !any(remaining)
            break
        end

        rem_idx = findall(remaining)
        num_rem = length(rem_idx)
        dW = sqrt(dt) * randn(rng, num_rem)
        X_rem = X[rem_idx]
        drift = r .* X_rem .* (1 .- X_rem / K)
        diff = sigma .* X_rem .* dW
        X_rem .+= drift .* dt .+ diff
        X_rem = max.(X_rem, 0.0)
        X[rem_idx] = X_rem

        # Harvest check
        harvested = X_rem .>= Wh
        if any(harvested)
            append!(harvest_times, fill(t, count(harvested)))
            remaining[rem_idx[harvested]] .= false
        end

        # Update rem_idx after harvest
        rem_idx = findall(remaining)
        num_rem = length(rem_idx)

        # Mortality
        if num_rem > 0
            deaths = rand(rng, num_rem) .< mu * dt
            if any(deaths)
                append!(death_times, fill(t, count(deaths)))
                remaining[rem_idx[deaths]] .= false
            end
        end

        # Track at specified times
        while length(results) < length(times) && t >= times[length(results) + 1]
            mean_weight = any(remaining) ? mean(X[remaining]) : NaN
            sd_weight = any(remaining) ? std(X[remaining]) : NaN
            H_mean_weight = any(harvested) ? mean(X_rem[harvested]) : 0.0
            H_sd_weight = any(harvested) ? std(X_rem[harvested]) : 0.0
            push!(results, (times[length(results) + 1], 
                            count(remaining), mean_weight, sd_weight,
                            count(harvested), H_mean_weight, H_sd_weight))
        end
    end

    # Fill remaining track times if early stop
    while length(results) < length(times)
        push!(results, (times[length(results) + 1], 0, NaN, NaN))
    end

    # println("Simulation results:")
    # for (time, rem, mean_w) in results
    #     mean_str = isnan(mean_w) ? "-" : string(round(mean_w, digits=2))
    #     println("Time: $time.0, Remaining: $rem, Mean Weight: $mean_str")
    # end

    # if !isempty(harvest_times)
    #     avg_harvest_time = mean(harvest_times)
    #     println("Average time to harvest (for harvested fish): $(round(avg_harvest_time, digits=2)) days")
    # else
    #     println("No fish harvested")
    # end

    # total_harvested = length(harvest_times)
    # total_died = length(death_times)
    # final_t = !any(remaining) ? times[end] : T  # Approximate
    # println("Total harvested: $total_harvested")
    # println("Total died: $total_died")
    # println("Tank empty at time: $(round(final_t, digits=1)) days")

    return results #, harvest_times, death_times
end

# Example call
tmp = simulate_popfish(50000, 0.01, 12.0, 0.05, 4.0, 0.001, 366.0, 7.0)
tmp = stack(tmp, dims=1)

p1 = plot(1:(T+1), MW, ribbon = 1.96 .* SW ./ sqrt(300), title = "Weekly harvested biomass")
p2 = plot(1:(T+1), cMW, ribbon = 1.96 .* cSW ./ sqrt(300), title = "Cumulative harvested biomass")
plot(p1, p2, layout=(1,2), size=(750,300))