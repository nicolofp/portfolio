
functions {
  real[] salmon_farm_ode(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    // State variables:
    // y[1]: number of salmon (N)
    // y[2]: mean weight (W)
    // y[3]: oxygen level (O)
    // y[4]: food per salmon (F)
    // y[5]: harvested salmon (H)
    
    // Parameters:
    // theta[1]: tank volume (m^3)
    // theta[2]: target weight (kg)
    // theta[3]: max food per salmon (kg/day)
    // theta[4]: food conversion ratio
    
    real N = y[1];  // Number of salmon
    real W = y[2];  // Mean weight (kg)
    real O = y[3];  // Oxygen level (mg/L)
    real F = y[4];  // Food per salmon (kg/day)
    real H = y[5];  // Harvested salmon (cumulative)
    
    real tank_volume = theta[1];
    real target_weight = theta[2];
    real max_food_per_salmon = theta[3];
    real food_conversion_ratio = theta[4];
    
    // Derived quantities
    real biomass_density = N * W / tank_volume;
    
    // Factors affecting growth
    real oxygen_factor = fmin(O / 8.0, 1.0);  // Normalize to optimal level (8 mg/L)
    real crowding_factor = exp(-0.05 * biomass_density);  // Exponential decrease with density
    real food_efficiency = F / (F + 0.01) * (1 - 0.1 * (W / target_weight));  // Efficiency decreases as fish grow
    
    // Growth rate calculation
    real growth_rate = 0.012 * oxygen_factor * crowding_factor * food_efficiency * F / food_conversion_ratio;
    
    // Oxygen dynamics
    real oxygen_consumption = 0.001 * biomass_density * (1 + 2 * F / max_food_per_salmon);
    real oxygen_replenishment = 0.1 * (10 - O);  // Approaches equilibrium of 10 mg/L
    
    // Food adjustment
    real food_adjustment = 0.01 * (max_food_per_salmon - F) * oxygen_factor;
    
    // Harvest rate - only harvest fish that reach target weight
    // Assuming 5% of fish above target weight are harvested daily
    real percent_above_target;
    if (W < target_weight) {
      percent_above_target = 0.05 * exp(-(target_weight - W)^2 / (2 * 0.5^2));
    } else {
      percent_above_target = 0.05;
    }
    real harvest_rate = N * percent_above_target;
    
    // Differential equations
    real dN_dt = -harvest_rate;  // Change in salmon number due to harvesting
    real dW_dt = growth_rate;    // Change in mean weight
    real dO_dt = oxygen_replenishment - oxygen_consumption;  // Change in oxygen
    real dF_dt = food_adjustment;  // Change in food per salmon
    real dH_dt = harvest_rate;   // Accumulated harvested salmon
    
    return {dN_dt, dW_dt, dO_dt, dF_dt, dH_dt};
  }
}

data {
  int<lower=1> T;                // Number of time points to solve for
  real t0;                       // Initial time point
  real ts[T];                    // Time points to solve for
  real<lower=0> y0[5];           // Initial state
  real<lower=0> theta[4];        // Parameters
}

transformed data {
  real x_r[0];                   // No real data for the ODE function
  int x_i[0];                    // No integer data for the ODE function
}

model {
  // Empty model block as we're just using Stan for the ODE solver
}

generated quantities {
  real y_hat[T, 5];              // Predicted states at each time point
  
  // Solve the ODE system
  y_hat = integrate_ode_rk45(
    salmon_farm_ode, y0, t0, ts, theta, x_r, x_i
  );
}

