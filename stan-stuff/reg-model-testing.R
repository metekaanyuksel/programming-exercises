setwd("~/coding-stuffs/stan-stuff")

require(rstan)
require(tidyverse)
require(posterior)
require(brms)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# https://mc-stan.org/rstan/articles/stanfit_objects.html
# https://mc-stan.org/rstanarm/index.html
# https://jrnold.github.io/bugs-examples-in-stan/index.html
# https://mc-stan.org/docs/functions-reference/index.html

### linear regression:

x <- runif(n = 500, min = 0, max = 1)
y <- rnorm(n = 500, mean = 2 + 3 * x, sd = 1)

data <- list(N = length(x), x = matrix(x, nrow = length(x)), y = y, K = 1)

# stan_model("linear-reg-model.stan") -> model
stan(file = "linear-reg-model.stan", data = data, iter = 2000) -> fit

summary(fit) -> summary_fit

# print(summary_fit)
# traceplot(fit)
# pairs(fit)

matrix_of_draws <- as.matrix(fit)

plot(matrix_of_draws[,1], matrix_of_draws[,2])
hist(matrix_of_draws[,1])

code <- get_stancode(fit)
cat(code)

inits <- get_inits(fit)
inits_chain1 <- inits[[1]]
print(inits_chain1)

print(get_elapsed_time(fit))

xgrid <- seq(0,1,length=100)

pdf("figures/bayes-lin-reg.pdf", width = 6, height = 4)

plot(data$x, data$y, pch = 1, col = "grey", xlab = "x", ylab = "y")

for (j in 1:50){
  
  alpha <- sample(x = matrix_of_draws[,1], size = 1)
  beta <- sample(x = matrix_of_draws[,2], size = 1)
  
  lines(xgrid, alpha + beta*xgrid, col = "black")
}

lines(xgrid, 
      as.numeric(summarise_draws(fit)[1,2]) + 
        as.numeric(summarise_draws(fit)[2,2])*xgrid, lwd = 5,
      col = "red")

dev.off()

### logistic regression with one predictor:

read_csv("~/coding-stuffs/data/farrell-davies-data.csv") -> disease_dataset
ifelse(disease_dataset$Deaths > 0, 1, 0) -> AnyDeaths

data <- list(y = AnyDeaths,
             x = disease_dataset$EvoIso,
             N = length(AnyDeaths)
)

# stan_model("logistic-reg-one-pred.stan") -> model_logistic_reg
stan(file = "logistic-reg-one-pred.stan", data = data, iter = 10000) -> fit_logistic

summary(fit_logistic)
pairs(fit_logistic)

summarise_draws(fit_logistic)

plot(fit_logistic, 
     show_density = TRUE, 
     ci_level = 0.95, 
     pars = "beta",
     fill_color = "grey")

hist(as.matrix(fit_logistic)[,2])

## plot 50 draws from posterior on top of the data

pdf("figures/farrell-davies-data.pdf", width = 6, height = 4)

plot(disease_dataset$EvoIso, AnyDeaths, 
     xlab = "evolutionary isolation", 
     ylab = "Pr(death)",
     pch = 19)

x <- seq(0,max(disease_dataset$EvoIso), length = 1000)

for (j in 1:50){
  
alpha <- sample(x = as.matrix(fit_logistic)[,1], size = 1)
beta <- sample(x = as.matrix(fit_logistic)[,2], size = 1)

lines(x,
  1/(1+exp(-alpha-beta*x)), 
   col = "lightgrey"
  )

}

alpha <- mean(as.matrix(fit_logistic)[,1])
beta <- mean(as.matrix(fit_logistic)[,2])

lines(x,
  1/(1+exp(-alpha-beta*x)), 
  col = "black", lwd = 3
)

# summarise_draws(fit_logistic)[1:2,"q5"] -> lower
# summarise_draws(fit_logistic)[1:2,"q95"] -> upper
# 
# lines(x,
#       1/(1+exp(-as.numeric(lower[1,])-as.numeric(lower[2,])*x)),
#       col = "black", lty = "dashed", lwd = 2
# )
# 
# lines(x,
#       1/(1+exp(-as.numeric(upper[1,])-as.numeric(upper[2,])*x)),
#       col = "black", lty = "dashed", lwd = 2
# )

dev.off()

### non-linear regression:

read_csv("~/coding-stuffs/data/chytrid-temp-data.csv") -> temperature_disease_data
temperature_disease_data$Zoospores <- round(temperature_disease_data$Zoospores)

temperature_disease_data |>
  ggplot(aes(x = Temperature, y = Zoospores, color = Study)) +
  geom_point(size = 3) +
  labs(x = "temperature", y = "# zoospores observed", color = "study") +
  theme_bw()

data <- list(N = dim(temperature_disease_data)[1],
             Y = temperature_disease_data$Zoospores,
             Temperature = temperature_disease_data$Temperature
             )

thermal_curve <- function(temp, Topt, w, A){
  return( A * exp( - ((temp - Topt) / (2 * w))^2 ) )
}

# 
# thermal_curve(temperature_disease_data$Temperature,
#               Topt = 18, w = 3, A = 5000) -> ZoosporeData
# 
# list(N = length(temperature_disease_data$Temperature),
#   Temperature = temperature_disease_data$Temperature, 
#   Y = rnbinom(n = length(temperature_disease_data$Temperature), 
#               mu = ZoosporeData, 
#               size = 0.3)
#   ) -> data
# 
# plot(data$Temperature, data$Y)

list() -> initial_conditions

for (i in 1:4){
  initial_conditions[[i]] <- list(A = 5000, 
                                  Topt = 18, 
                                  w = 3, 
                                  shape = 0.3)
}

stan(file = "nonlinear-reg.stan", data = data,
     warmup = 1000,
     iter = 10000,
     init = initial_conditions) -> fit

summarise_draws(fit)
pairs(fit)

plot(fit, 
     show_density = TRUE, 
     ci_level = 0.95, 
     # pars = c("Topt", "w", "shape", "A"),
     fill_color = "grey")

hist(as.matrix(fit)[,1])
plot(as.matrix(fit)[,1], as.matrix(fit)[,2])

write_csv(as.data.frame( as.matrix(fit) ), "posterior-draws-TPC.csv")

pdf("figures/TPC-chytrid.pdf", width = 6, height = 4)

plot(data$Temperature, data$Y, pch = 19, 
     xlab = "temperature", ylab = "# zoospores produced")

temp <- seq(0,max(data$Temperature), length = 1000)

for (j in 1:100){
  
  Topt <- sample(x = as.matrix(fit)[,1], size = 1)
  w <- sample(x = as.matrix(fit)[,2], size = 1)
  shape <- sample(x = as.matrix(fit)[,3], size = 1)
  A <- sample(x = as.matrix(fit)[,4], size = 1)
  
  lines(temp,
        thermal_curve(temp, Topt = Topt, w = w, A = A), 
        col = "lightgrey"
  )
  
}

points(data$Temperature, data$Y, pch = 19)

as.numeric(summarise_draws(fit)[1,2]) -> Topt_mean
as.numeric(summarise_draws(fit)[2,2]) -> w_mean
as.numeric(summarise_draws(fit)[3,2]) -> shape_mean
as.numeric(summarise_draws(fit)[4,2]) -> A_mean

lines(temp,
      thermal_curve(temp, Topt = Topt_mean, w = w_mean, A = A_mean), 
      col = "black",
      lwd = 3
)

dev.off()