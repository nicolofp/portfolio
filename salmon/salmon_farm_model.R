# Salmon Farm Model using cmdstanr
library(cmdstanr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(posterior)

stan_file = "salmon_final.stan"

# Function to initialize and simulate the salmon farm
simulate_salmon_farm <- function(
    initial_salmon = 50000,         # Initial number of salmon
    initial_mean_weight = 4.5,      # Initial mean weight (kg)
    tank_volume = 1000000,          # Tank volume (m^3)
    initial_oxygen = 2.0,           # Initial dissolved oxygen (mg/L)
    target_weight = 7.0,            # Target harvest weight (kg)
    max_food_per_salmon = 0.89,     # Maximum food per salmon per day (kg)
    food_conversion_ratio = 1.9,    # Food conversion ratio (food to weight gain)
    days = 365                      # Simulation duration
) {
  # Compile the Stan model
  mod <- cmdstan_model(stan_file)
  
  # Initial food amount per salmon (kg/day)
  initial_food_per_salmon <- max_food_per_salmon * 0.8  # Start at 80% of maximum
  
  # Initial state
  y0 <- c(
    initial_salmon,            # Number of salmon
    initial_mean_weight,       # Mean weight
    initial_oxygen,            # Oxygen level
    initial_food_per_salmon,   # Food per salmon
    0                          # Harvested salmon (cumulative)
  )
  
  # Parameters
  theta <- c(
    tank_volume,               # Tank volume
    target_weight,             # Target weight
    max_food_per_salmon,       # Max food per salmon
    food_conversion_ratio      # Food conversion ratio
  )
  
  # Time points for simulation
  ts <- seq(from = 1, to = days, by = 1)
  
  # Prepare data for Stan
  stan_data <- list(
    T = length(ts),
    ts = ts,
    y0 = y0,
    tank_volume = tank_volume,     
    target_weight = target_weight,
    max_food_per_salmon = max_food_per_salmon, 
    food_conversion_ratio = food_conversion_ratio 
  )
  
  # Run the simulation using fixed_param method
  fit <- mod$sample(
    data = stan_data,
    chains = 1,
    iter_sampling = 1,
    iter_warmup = 0,
    fixed_param = TRUE,
    refresh = 0
  )
  
  # Extract results
  draws_df <- fit$draws(format = "draws_df")
  
  # Process the results 
  # In cmdstanr, the array parameters are named y_hat[t,i] instead of y_hat[t,i,j]
  farm_data <- data.frame(
    time = ts,
    num_salmon = numeric(length(ts)),
    mean_weight = numeric(length(ts)),
    oxygen = numeric(length(ts)),
    food_per_salmon = numeric(length(ts)),
    harvested = numeric(length(ts))
  )
  
  for (t in 1:length(ts)) {
    farm_data$num_salmon[t] <- draws_df[[paste0("y_hat[", t, ",1]")]]
    farm_data$mean_weight[t] <- draws_df[[paste0("y_hat[", t, ",2]")]]
    farm_data$oxygen[t] <- draws_df[[paste0("y_hat[", t, ",3]")]]
    farm_data$food_per_salmon[t] <- draws_df[[paste0("y_hat[", t, ",4]")]]
    farm_data$harvested[t] <- draws_df[[paste0("y_hat[", t, ",5]")]]
  }
  
  # Add derived quantity - biomass
  farm_data <- farm_data %>%
    mutate(biomass = num_salmon * mean_weight)
  
  return(farm_data)
}

# Generate initial weights for visualization
simulate_initial_weights <- function(n = 1000, mean_weight = 2.5, sd_weight = 0.3) {
  weights <- rnorm(n, mean = mean_weight, sd = sd_weight)
  weights <- pmax(0.5, pmin(weights, mean_weight * 2))  # Clip to sensible range
  
  ggplot(data.frame(weight = weights), aes(x = weight)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    geom_vline(xintercept = mean_weight, linetype = "dashed", color = "red") +
    labs(title = "Initial Weight Distribution", x = "Weight (kg)", y = "Count") +
    theme_minimal()
}

# Function to plot results
plot_farm_results <- function(data) {
  p1 <- ggplot(data, aes(x = time, y = num_salmon)) +
    geom_line() +
    #geom_point() +
    labs(title = "Number of Salmon", x = "Days", y = "Count") +
    theme_minimal()
  
  p2 <- ggplot(data, aes(x = time, y = mean_weight)) +
    geom_line() +
    #geom_point() +
    geom_hline(yintercept = 7, linetype = "dashed", color = "red") +
    labs(title = "Mean Salmon Weight", x = "Days", y = "Weight (kg)") +
    theme_minimal()
  
  p3 <- ggplot(data, aes(x = time, y = oxygen)) +
    geom_line() +
    #geom_point() +
    labs(title = "Oxygen Level", x = "Days", y = "Dissolved Oxygen (mg/L)") +
    theme_minimal()
  
  p4 <- ggplot(data, aes(x = time, y = food_per_salmon)) +
    geom_line() +
    #geom_point() +
    labs(title = "Food per Salmon", x = "Days", y = "kg/day") +
    theme_minimal()
  
  p5 <- ggplot(data, aes(x = time, y = biomass)) +
    geom_line() +
    #geom_point() +
    labs(title = "Total Biomass", x = "Days", y = "kg") +
    theme_minimal()
  
  p6 <- ggplot(data, aes(x = time, y = harvested)) +
    geom_line() +
    #geom_point() +
    labs(title = "Cumulative Harvested Salmon", x = "Days", y = "Count") +
    theme_minimal()
  
  grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 2)
}

# Run the simulation (this takes a few seconds)
set.seed(42)  # For reproducibility
sim_results <- simulate_salmon_farm(days = 500)

# Display the simulation results
plot_farm_results(sim_results)

# Display initial weight distribution
weight_dist <- simulate_initial_weights()
print(weight_dist)
