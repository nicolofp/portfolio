functions {
  vector salmon_farm(real t,
                     vector y,
                     real tank_volume,    
                     real target_weight,
                     real max_food_per_salmon,
                     real food_conversion_ratio) { 
                    
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
    vector[5] dydt;
    
    // Derived quantities
    real biomass_density = y[1] * y[2] / tank_volume;
    
    // Factors affecting growth
    real oxygen_factor = fmin(y[3] / 8.0, 1.0);  // Normalize to optimal level (8 mg/L)
    real crowding_factor = exp(-0.05 * biomass_density);  // Exponential decrease with density
    real food_efficiency = y[4] / (y[4] + 0.01) * (1 - 0.1 * (y[2] / target_weight));  // Efficiency decreases as fish grow
    
    // Growth rate calculation
    real growth_rate = 0.012 * oxygen_factor * crowding_factor * food_efficiency * y[4] / food_conversion_ratio;
    
    // Oxygen dynamics
    real oxygen_consumption = 0.001 * biomass_density * (1 + 2 * y[5] / max_food_per_salmon);
    real oxygen_replenishment = 0.1 * (10 - y[3]);  // Approaches equilibrium of 10 mg/L
    
    // Food adjustment
    real food_adjustment = 0.01 * (max_food_per_salmon - y[4]) * oxygen_factor;
    
    // Harvest rate - only harvest fish that reach target weight
    real percent_above_target;
    if (y[2] < target_weight) {
      percent_above_target = 0.05 * exp(-(target_weight - y[2])^2 / (2 * 0.5^2));
    } else {
      percent_above_target = 0.05;
    }
    real harvest_rate = y[1] * percent_above_target;

    dydt[1] = -harvest_rate;  // Change in salmon number due to harvesting
    dydt[2] = growth_rate;    // Change in mean weight
    dydt[3] = oxygen_replenishment - oxygen_consumption;  // Change in oxygen
    dydt[4] = food_adjustment;  // Change in food per salmon
    dydt[5] = harvest_rate;   // Accumulated harvested salmon

    return dydt;
  }
}

data {
  int<lower=1> T;                // Number of time points to solve for
  array[T] real ts;              // Time points to solve for
  vector[5] y0;
  real tank_volume;    
  real target_weight;
  real max_food_per_salmon;
  real food_conversion_ratio;    
}

transformed data {
  real t0 = 0; 
}

model {
  // Empty model block as we're just using Stan for the ODE solver
}

generated quantities {
  // Solve the ODE system
  array[T] vector[5] y_hat = ode_bdf(salmon_farm, y0, t0, ts, tank_volume, target_weight, max_food_per_salmon, food_conversion_ratio);
  //for (t in 1:T) {
  //  y_hat[t, 1] += normal_rng(0, 0.1);
  //  y_hat[t, 2] += normal_rng(0, 0.1);
  //  y_hat[t, 3] += normal_rng(0, 0.1);
  //  y_hat[t, 4] += normal_rng(0, 0.1);
  //  y_hat[t, 5] += normal_rng(0, 0.1);
  //}
}
