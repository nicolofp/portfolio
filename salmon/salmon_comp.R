# Salmon Farm Model using cmdstanr
library(cmdstanr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(posterior)

stan_file = "salmon_new.stan"
mod <- cmdstan_model(stan_file)

# Initial state
y0 <- c(10000,
        3,
        0.5,
        0,
        0,
        0.1)

# Time points for simulation
ts <- seq(from = 1, to = 52, by = 1)

# Prepare data for Stan
stan_data <- list(
  T = length(ts),
  ts = ts,
  y0 = y0,
  h = 1, 
  gr = 0.05, 
  tw = 5,
  mu = 0.001, 
  W_max = 12, 
  sigma_g = 0.1
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
farm_data <- data.frame(
  time = ts,
  num_salmon = numeric(length(ts)),
  mean_weight = numeric(length(ts)),
  sd_weight = numeric(length(ts)),
  har_salmon = numeric(length(ts)),
  mean_weight_har_salmon = numeric(length(ts)),
  sd_weight_har_salmon = numeric(length(ts))
)

for (t in 1:length(ts)) {
  farm_data$num_salmon[t] <- draws_df[[paste0("y_hat[", t, ",1]")]]
  farm_data$mean_weight[t] <- draws_df[[paste0("y_hat[", t, ",2]")]]
  farm_data$sd_weight[t] <- draws_df[[paste0("y_hat[", t, ",3]")]]
  farm_data$har_salmon[t] <- draws_df[[paste0("y_hat[", t, ",4]")]]
  farm_data$mean_weight_har_salmon[t] <- draws_df[[paste0("y_hat[", t, ",5]")]]
  farm_data$sd_weight_har_salmon[t] <- draws_df[[paste0("y_hat[", t, ",6]")]]
}

# Add derived quantity - biomass
farm_data <- farm_data %>%
  mutate(biomass = num_salmon * mean_weight,
         harvested = har_salmon * mean_weight_har_salmon)


ggplot(farm_data, aes(x = time, y = harvested)) +
  geom_line() +
  geom_line(aes(x = time, y = biomass)) +
  #geom_point() +
  labs(title = "Biomass over time", x = "Days", y = "Biomass") +
  theme_minimal()

ggplot(farm_data, aes(x = time, y = mean_weight_har_salmon*har_salmon )) +
  geom_line() +
  geom_ribbon(aes(ymin = (har_salmon*mean_weight_har_salmon  - 1.96*har_salmon*sd_weight_har_salmon), 
                  ymax = (har_salmon*mean_weight_har_salmon  + 1.96*har_salmon*sd_weight_har_salmon)), 
              fill = "grey70", alpha = 0.4) + 
  #geom_point() +
  labs(title = "Biomass over time", x = "Days", y = "Biomass") +
  theme_minimal()

ggplot(farm_data, aes(x = time, y = mean_weight_har_salmon)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_weight_har_salmon  - 1.96*sd_weight_har_salmon, 
                  ymax = mean_weight_har_salmon  + 1.96*sd_weight_har_salmon), 
              fill = "grey70", alpha = 0.4) + 
  #geom_point() +
  labs(title = "Biomass over time", x = "Days", y = "Biomass") +
  theme_minimal()

ggplot(farm_data, aes(x = time, y = mean_weight )) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_weight   - 1.96*sd_weight , 
                  ymax = mean_weight   + 1.96*sd_weight ), 
              fill = "grey70", alpha = 0.4) + 
  #geom_point() +
  labs(title = "Biomass over time", x = "Days", y = "Biomass") +
  theme_minimal()


calculate_W_inst <- function(W, V) {
  # W: Mean weight
  # V: Variance of weight
  # Returns W_inst: Instantaneous weight
  
  alpha <- (5 - W) / sqrt(V)
  
  # Calculate phi(alpha) - PDF of standard normal at alpha
  phi_alpha <- dnorm(alpha)
  
  # Calculate 1-Phi(alpha) - Complement of CDF of standard normal at alpha
  complement_Phi_alpha <- 1 - pnorm(alpha)
  
  # Calculate W_inst using the formula
  W_inst <- W + sqrt(V) * (phi_alpha / complement_Phi_alpha)
  
  return(W_inst)
}

calculate_V_inst <- function(W, V) {
  # W: Mean weight
  # V: Variance of weight
  # Returns V_inst: Instantaneous variance
  
  alpha <- (5 - W) / sqrt(V)
  
  # Calculate phi(alpha) - PDF of standard normal at alpha
  phi_alpha <- dnorm(alpha)
  
  # Calculate 1-Phi(alpha) - Complement of CDF of standard normal at alpha
  complement_Phi_alpha <- 1 - pnorm(alpha)
  
  # Calculate the ratio term
  ratio <- phi_alpha / complement_Phi_alpha
  
  # Calculate V_inst using the formula
  term1 <- (phi_alpha * alpha) / complement_Phi_alpha
  term2 <- (ratio)^2
  V_inst <- V * (1 + term1 - term2)
  
  return(V_inst)
}


# growth weekly 2.5%
set.seed(1990)
x = vector("list",50)
x[[1]] = rnorm(10000,2,0.5)
rs = data.table(N = c(NROW(x[[1]]),rep(NA,49)),
                M = c(mean(x[[1]]),rep(NA,49)),
                S = c(sd(x[[1]]),rep(NA,49)))

for(i in 2:50){
  tmp = x[[i-1]]*(1 + 0.025)
  x[[i]] = tmp[tmp < 5]
  rs$N[i] <- NROW(tmp)
  rs$M[i] <- mean(tmp)
  rs$S[i] <- sd(tmp)
}

plot(1:50,rs$M,type="l")
plot(1:50,rs$S,type="l")

