library(deSolve)

# =============================================================================
# Parameters (adjust these as needed)
# =============================================================================
n <- 200                    # Number of weight bins
w_min <- 0.01               # Minimum weight
w_max <- 5.5                # Maximum weight
a <- 0.5                    # Growth parameter: g(w) = a - b * w
b <- 0.1                    # Growth parameter
w_harvest <- 3.5            # Harvest threshold weight
H <- 0.8                    # Harvesting rate for w >= w_harvest
mu <- 0.05                  # Natural mortality rate (constant)
N0 <- 50000                 # Initial total population
meanlog <- log(0.1)         # Lognormal meanlog for initial distribution
sdlog <- 0.4                # Lognormal sdlog
t_max <- 52                 # Simulation time
n_times <- 52*3.5           # Number of time points

# Weight vector (bin centers)
w <- seq(w_min, w_max, length.out = n)
dw <- w[2] - w[1]

# Growth rate vector: g(w)
g_vec <- a - b * w

# Mortality and harvesting vectors
mu_vec <- rep(mu, n)
h_vec <- ifelse(w >= w_harvest, H, 0)

# Initial condition: normalized lognormal density + Bh=0
u0_cont <- dlnorm(w, meanlog = meanlog, sdlog = sdlog)
integral_u0 <- sum(u0_cont * dw)  # Approximate integral
u0 <- c((u0_cont / integral_u0) * N0, 0)  # Append Bh=0

# ODE derivative function (method of lines + harvested biomass)
model <- function(t, u, parms) {
  g_vec <- parms$g_vec
  mu_vec <- parms$mu_vec
  h_vec <- parms$h_vec
  dw <- parms$dw
  n <- parms$n
  w <- parms$w
  
  du <- numeric(n + 1)  # Extra for Bh
  # u bins
  du[1] <- -g_vec[1] * u[1] / dw - (mu_vec[1] + h_vec[1]) * u[1]
  for (i in 2:n) {
    du[i] <- (g_vec[i-1] * u[i-1] - g_vec[i] * u[i]) / dw - (mu_vec[i] + h_vec[i]) * u[i]
  }
  # Harvested biomass increment: sum h_i * u_i * w_i * dw
  du[n+1] <- sum(h_vec * u[1:n] * w * dw)
  list(du)
}

# Parameters list for ode()
parms <- list(g_vec = g_vec, mu_vec = mu_vec, h_vec = h_vec, dw = dw, n = n, w = w)

# Time points
times <- seq(0, t_max, length.out = n_times)

# Solve the ODE system
out <- ode(y = u0, times = times, func = model, parms = parms)

# =============================================================================
# Outputs and Plots
# =============================================================================
# Extract u(w,t) and Bh(t)
u_out <- out[, 2:(n+1)]  # Densities
Bh <- out[, n+2]         # Cumulative harvested biomass

# Total population and biomass over time (approximate integrals)
totalN <- rowSums(u_out * dw)
totalB <- rowSums(sweep(u_out, 2, w, "*") * dw)  # sum w u dw

# Incremental harvested biomass between timepoints
delta_Bh <- c(0, diff(Bh))  # Harvested in each interval, starting with 0 at t=0

# Plots
par(mfrow = c(2, 2))
# Total population
plot(times, totalN, type = "l", lwd = 2, main = "Total Population", xlab = "Time", ylab = "N(t)", col = "blue")
# Total biomass
plot(times, totalB, type = "l", lwd = 2, main = "Total Biomass", xlab = "Time", ylab = "B(t)", col = "green")
# Cumulative harvested biomass
plot(times, Bh, type = "l", lwd = 2, main = "Cumulative Harvested Biomass", xlab = "Time", ylab = "Bh(t)", col = "red")
# Incremental harvested biomass
plot(times, delta_Bh, type = "l", lwd = 2, main = "Harvested Biomass per Interval", xlab = "Time", ylab = "ΔBh(t)", col = "orange")
par(mfrow = c(1, 1))

# Optional: Print harvested at each timepoint
harvested_df <- data.frame(time = times, cumulative_Bh = Bh, incremental_Bh = delta_Bh)
print(head(harvested_df))  # Show first few for example


# Density distributions: initial vs final
par(mfrow = c(1, 2))
plot(w, u_out[1,], type = "l", lwd = 2, 
     main = "Initial u(w)", xlab = "Weight w", ylab = "u(w)")
plot(w, u_out[nrow(u_out),], type = "l", lwd = 2, 
     main = "Final u(w)", xlab = "Weight w", ylab = "u(w)")
par(mfrow = c(1, 1))

# Optional: Heatmap of u(w,t) evolution
image(x = times, y = w, z = u_out, 
      xlab = "Time", ylab = "Weight w", main = "u(w,t) Evolution",
      col = heat.colors(100))


plot(w, u_out[1, ], type = "l", #lwd = 2, 
     main = "Initial u(w)", xlab = "Weight w", ylab = "u(w)")
for(i in 2:nrow(u_out)){
  lines(w, u_out[i,], type = "l")
}