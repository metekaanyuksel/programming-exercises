require(tidyverse)
require(rstan)
require(posterior)
require(deSolve)

# https://blog.djnavarro.net/posts/2023-05-16_stan-ode/#when-biology-isnt-analytically-tractable
# https://www.magesblog.com/post/2021-02-08-fitting-multivariate-ode-models-with-brms/
# https://peter-stewart.github.io/blog/classic-ecological-models-exponential-logistic-growth/

setwd("~/coding-stuffs/stan-stuff")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

LV_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters,time)), {
    
   du = (alpha - beta * v) * u
   dv = (-gamma + delta * u) * v
    
    return(list(c(du, dv)))
  })
}

# params <- c(alpha = 0.545, gamma = 0.028, beta = 0.803, delta = 0.024)
# state <- c(u=33.956, v=5.933)
# times <- seq(0, 1000, 0.1)
# out <- as.data.frame(ode(state, times, LV_model, params))
# plot(v~time, out, type = "l")

### inference of Lotka-Volterra parameters using hare-lynx data:

read.csv("~/coding-stuffs/data/hudson-bay-lynx-hare.csv",
         comment.char="#") -> dat

N <- length(dat$year) - 1
ts <- 1:N
y_init <- c(dat$hare[1], dat$lynx[1])
y <- as.matrix(dat[2:(N + 1), 2:3])
y <- cbind(y[ , 2], y[ , 1]); # hare, lynx order

lynx_hare_data <- list(N = N, ts = ts, y_init = y_init, y = y)

model <- stan_model("Lotka-Volterra-fitting.stan")
fit <- sampling(model, data = lynx_hare_data, seed = 123)

summarise_draws(fit)

write_csv(as.data.frame( as.matrix(fit) ), "posterior-draws-LV.csv")

pdf("figures/LV-fit.pdf", width = 6, height = 4)

plot(lynx_hare_data$ts, lynx_hare_data$y[,1], type = "l",
     xlab = "time", ylab = "# hares in black, lynx in red", ylim = 
       c(0,max(max(lynx_hare_data$y[,1]), max(lynx_hare_data$y[,2])))
     )

lines(lynx_hare_data$ts, lynx_hare_data$y[,2], type = "l",
     xlab = "time", ylab = "# lynx", col = "red")

times <- seq(0,max(lynx_hare_data$ts), length = 1000)

for (j in 1:20){
  
  alpha <- sample(x = as.matrix(fit)[,1], size = 1)
  beta <- sample(x = as.matrix(fit)[,2], size = 1)
  gamma <- sample(x = as.matrix(fit)[,3], size = 1)
  delta <- sample(x = as.matrix(fit)[,4], size = 1)
  initialH <- sample(x = as.matrix(fit)[,5], size = 1)
  initialL <- sample(x = as.matrix(fit)[,6], size = 1)
  
  params <- c(alpha = alpha, 
              gamma = gamma, 
              beta = beta, 
              delta = delta)
  
  state <- c(u=initialH, v=initialL)
  
  out <- as.data.frame(ode(state, times, LV_model, params))
  
  lines(out$time,
        out$u, 
        col = "lightgrey"
  )
  
  lines(out$time,
        out$v, 
        col = "pink"
  )
  
}

lines(lynx_hare_data$ts, lynx_hare_data$y[,1], type = "l",
     xlab = "time", ylab = "# hares in black, lynx in red", ylim = 
       c(0,max(max(lynx_hare_data$y[,1]), max(lynx_hare_data$y[,2])))
)

lines(lynx_hare_data$ts, lynx_hare_data$y[,2], type = "l",
      xlab = "time", ylab = "# lynx", col = "red")

as.numeric(summarise_draws(fit)[1,2]) -> alpha
as.numeric(summarise_draws(fit)[2,2]) -> beta
as.numeric(summarise_draws(fit)[3,2]) -> gamma
as.numeric(summarise_draws(fit)[4,2]) -> delta
as.numeric(summarise_draws(fit)[5,2]) -> initialH
as.numeric(summarise_draws(fit)[6,2]) -> initialL

params <- c(alpha = alpha, 
            gamma = gamma, 
            beta = beta, 
            delta = delta)

state <- c(u=initialH, v=initialL)

out <- as.data.frame(ode(state, times, LV_model, params))

lines(out$time,
      out$u, 
      col = "black",
      lwd = 3
)

lines(out$time,
      out$v, 
      col = "red",
      lwd = 3
)

dev.off()