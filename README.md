# gompertz-jags
This implements the Gompertz distribution in JAGS. It follows the instructions published by D. Wabersich and J. Vandekerckhove (2013). DOI: https://doi.org/10.3758/s13428-013-0369-3

To install, issue the following commands in your terminal:

```
autoreconf -fvi
./configure
make
sudo make install
```

To use it in JAGS, type `load.module("Gompertz")` and in your model file e. g. `dgomp(0.05, 0.005)`.

Here is a workable example:

```
library(rjags)
library(mortAAR)
pop_sim <- mortAAR::pop.sim.gomp(100,M = 50, start_age = 15)

modelString = "
model {
  for (i in 1:N) {
    y[i] ~ dgomp(b, a)
  }
  b ~ dunif(0.02, 0.1)
  a ~ dnorm(0,1) T(0,0.02)
}
" # close quote for modelString

# Path to the custom JAGS module
load.module("Gompertz")

# Data
data <- list(N = nrow(pop_sim$result), y = pop_sim$result$age - 15)

# JAGS model file
model_file <- textConnection(modelString)

# Parameters to monitor
parameters <- c("b", "a")

# Initialize the JAGS model
model <- jags.model( model_file, data = data, n.chains = 3, n.adapt = 1000)

# Burn-in period
update(model, n.iter = 2000)

# Sample from the posterior
samples <- coda.samples(model, variable.names = parameters, n.iter = 10000)

sampled <- samples[[1]]
print(summary(sampled))
plot(samples)
```
