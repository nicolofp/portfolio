functions {
  vector ode_system(real t, 
                    vector y,
                    real h,
                    real gr, 
                    real mu,
                    real tw,
                    real W_max,
                    real sigma_g) {
    // State variables:
    // y[1]: total number of salmon at time t
    // y[2]: mean individual salmon weight (kg) at time t
    // y[3]: variance of salmon weight (kg2), representing the spread of the size distribution.
    
    // Parameters:
    // h      : % of harvested fish over treshold
    // gr     : growth ratio (%)
    // tw     : target weight (kg)
    // mu     : death rate (%)
    // W_max  : max weight allowed per fish
    // sigma_g: sd of growth rate  
    
    vector[6] dydt;
    
    // Harvesting calculations
    real alpha = y[3] > 0 ? (tw - y[2]) / sqrt(y[3]) : 0;
    real phi_alpha = exp(-0.5 * square(alpha)) / sqrt(2 * pi());
    real Phi_alpha = y[3] > 0 ? Phi(alpha) : 1;
    real harvest_frac = y[3] > 0 ? 1 - Phi_alpha : 0;
    
    // Post-harvest mean and variance
    real W_post = y[2] - (y[3] > 0 ? sqrt(y[3]) * phi_alpha / fmax(Phi_alpha, 1e-6) : 0);
    real V_post = y[3] * (1 - (y[3] > 0 ? phi_alpha * alpha / fmax(Phi_alpha, 1e-6) + square(phi_alpha / fmax(Phi_alpha, 1e-6)) : 0));
    
    // Instantaneous harvested mean and variance
    real W_inst = y[2] + (y[3] > 0 ? sqrt(y[3]) * phi_alpha / fmax(1 - Phi_alpha, 1e-6) : 0);
    real V_inst = y[3] * (1 + (y[3] > 0 ? phi_alpha * alpha / fmax(1 - Phi_alpha, 1e-6) - square(phi_alpha / fmax(1 - Phi_alpha, 1e-6)) : 0));
    
    // ODEs
    dydt[1] = -mu * y[1] - h * harvest_frac * y[1];
    dydt[2] = gr * (1 - y[2] / W_max) - h * harvest_frac * (y[2] - W_post);
    dydt[3] = square(sigma_g) * gr - h * harvest_frac * (y[3] - V_post);
    dydt[4] = h * harvest_frac * y[1];
    dydt[5] = h * harvest_frac * y[1] * (W_inst - y[5]) / fmax(y[4] + 1e-6, 1e-6); 
    dydt[6] = h * harvest_frac * y[1] * (V_inst + square(W_inst - y[5]) - y[6]) / fmax(y[4] + 1e-6, 1e-6);
    
    return dydt;
  }
}

data {
  int<lower=1> T;                // Number of time points to solve for
  array[T] real ts;              // Time points to solve for
  vector[6] y0;
  real h;    
  real gr;
  real tw;
  real mu;
  real W_max;
  real<lower=0> sigma_g;
}

transformed data {
  real t0 = 0; 
}

model {
  // Empty model block as we're just using Stan for the ODE solver
}

generated quantities {
  // Solve the ODE system
  array[T] vector[6] y_hat = ode_bdf(ode_system, y0, t0, ts, h, gr, mu, tw, W_max, sigma_g);
  // for (t in 1:T) {
  //   y_hat[t, 1] += normal_rng(0, 0.5);
  //   y_hat[t, 2] += normal_rng(0, 0.5);
  //   y_hat[t, 3] += normal_rng(0, 0.5);
  //   y_hat[t, 4] += normal_rng(0, 0.5);
  //   y_hat[t, 5] += normal_rng(0, 0.5);
  // }
}
