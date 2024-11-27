# gompertz-jags
This implements the Gompertz distribution in JAGS. It follows the instructions published by D. Wabersich and J. Vandekerckhove (2013). DOI: https://doi.org/10.3758/s13428-013-0369-3

Beside the scale and shape parameter, it includes a parameter for maximum age because otherwise unrealistic high ages (at least for humans) will be sampled.
To install, issue the following commands in your terminal:

autoreconf -fvi

./configure

make

sudo make install

To use it in JAGS, type `load.module("Gompertz")` and in your model file e. g. `dgomp(0.05, 0.005, 85)` for a maximum age of 85 (remember that the Gompertz function by definition starts at the relative age of 0. With a real starting age of 15 (as is often used), the maximum age would be 100).
