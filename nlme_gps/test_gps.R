library(data.table)
library(bayesplot)
library(cmdstanr)
library(jsonlite)
library(ggplot2)
library(brms)

list_data = list.files("G:/My Drive/Temporary/gps_20250419/activities/",full.names = T)

data_training = lapply(list_data, function(i){
  print(i)
  x = read_json(i)
  if(x$latitude == 0 && !is.null(x$latitude) && !is.null(x$datastream)){
    tmp = rbindlist(x$datastream, fill = TRUE, 
                    use.names = TRUE, idcol = "id")
    tmp$date = as.Date(substr(x$startdate,1,10))  
    tmp = tmp[,.(date,id = as.numeric(id),hr = `1`)]
  }else{ NULL }
})
data_training = rbindlist(data_training)
data_training_red = data_training[date > "2022-12-04" & id > 1780]
data_training_red[, hr_max := max(hr), by = date]
data_training_red[, is_hr_max := ifelse(hr == hr_max,1,0)]

training_date = data_training_red[,.(which(is_hr_max == 1),min(hr)), 
                              by = date][V1 < 10 & date != "2023-03-18" & date != "2023-03-20", 
                                         date]

dt = data_training_red[date %in% training_date]
dt = dt[date %in% dt[id == 1785 & hr > 150,date]]
dt = dt[!(date %in% dt[id < 2250 & hr < 75,date])]

# Covariates month and number of running session in the previous 2w
cov1 = lapply(sort(unique(dt$date)), function(i){
  tmp = sum(sort(unique(dt$date)) %between% c(as.Date(i-14),as.Date(i-1)))
  data.table(date = i, month = month(i), N_train = tmp)
})
cov1 = rbindlist(cov1)
dt = merge(dt,cov1,by = "date")
dt = dt[id >= 1800]
dt[, time := (id-1800)/60]
dt[,':='(month = as.factor(month))]

dt_test = dt[time %in% seq(0,15,0.5)]  
date_exclude = dt_test[,.N, by = date][N<10,date]
dt_test = dt_test[!(date %in% date_exclude)]


ggplot(dt_test) +
  geom_line(aes(x = time, y = hr, group = date), alpha = 0.3) +
  geom_function(fun = function(x) 95 * exp(-0.32*x) + 92.43, 
                col = "red", linewidth = 1)

fit_exponential <- brms::brm(
  formula = brms::bf(hr ~ a * exp(-exp(logb) * time) + k, 
                     a + logb ~ 1 + N_train + (1 | date), 
                     k ~ 1,
                     nl=TRUE),
  data = dt_test,
  cores = 4,
  backend = "cmdstan",
  init = replicate(4,list(b_a_Intercept = 95       # Fixed intercept for a
                          ,b_a_N_train = 0         # Effect of N_train on a
                          ,b_b_Intercept = 0
                          ,b_b_N_train = 0
                          ,b_k_Intercept = 90
                          ), 
                   simplify = F),
  seed = 85538,
  prior   = c(prior(normal(96,1), nlpar = "a", coef="Intercept"),
              prior(normal(0,1), nlpar = "a", coef="N_train"),
              prior(normal(0,1), nlpar = "logb", coef="Intercept"),
              prior(normal(0,1), nlpar = "logb", coef="N_train"),
              prior(normal(90,1), nlpar = "k", lb = 0)),
  control = list(adapt_delta = 0.99,
                 max_treedepth = 15),
  file = "fit_exp_re_cov.rds"
)

mcmc_acf(fit_exponential)
mcmc_dens(fit_exponential)
mcmc_trace(fit_exponential)

bayes_R2(fit_exponential)
bayes_R2(fit_power)


saveRDS(fit_exponential,"fit_exp_good1.rds")
plot(conditional_effects(fit_exponential), ask = FALSE)
mcmc_plot(fit_exponential)


dt_test = dt[date %in% unique(dt$date)[c(6)] & 
  time %in% seq(0,15,0.5)]  
Is = ranef(fit_exponential)$date[6,1,1]
Ks = exp(-1.01 - 0.01 + ranef(fit_exponential)$date[6,1,2])

ggplot(dt_test) +
  geom_line(aes(x = time, y = hr, group = date), alpha = 0.3) +
  geom_function(fun = function(x) 93.16 * exp(-0.3828*x) + 92.74, 
                col = "red", linewidth = 1) +
  # geom_function(fun = function(x) (94.99 + Is) * exp(-(Ks)*x) + 92.43, 
  #               col = "blue") + 
  geom_function(fun = function(x) 10.55 * x^(-0.236) + 87.86, 
                col = "blue", linewidth = 1)

# a≈10.55, b≈0.236, c≈87.86
# https://michael-franke.github.io/Bayesian-Regression/practice-sheets/10a-nonLinear.html
# https://cran.r-project.org/web/packages/brms/vignettes/brms_nonlinear.html
# https://cran.r-project.org/web/packages/brms/vignettes/brms_distreg.html


fit_power <- brms::brm(
  formula = brms::bf(hr ~ a * time^exp(-logb) + k, 
                     a + logb ~ 1 + N_train + (1 | date), 
                     k ~ 1,
                     nl=TRUE),
  data = dt_test,
  cores = 4,
  backend = "cmdstan",
  init = replicate(4,list(b_a_Intercept = 10       # Fixed intercept for a
                          ,b_a_N_train = 0         # Effect of N_train on a
                          ,b_logb_Intercept = 1
                          ,b_logb_N_train = 0
                          ,b_k_Intercept = 90
                          ), simplify = F),
  seed = 1243547,
  prior   = c(prior(normal(10,1), nlpar = "a", coef="Intercept"),
              prior(normal(0,1), nlpar = "a", coef="N_train"),
              prior(normal(0,1), nlpar = "logb", coef="Intercept"),
              prior(normal(0,1), nlpar = "logb", coef="N_train"),
              prior(normal(90,1), nlpar = "k", lb = 0)),
  control = list(adapt_delta = 0.99,
                 max_treedepth = 15),
  file = "fit_pow_re_cov_test.rds"
)

