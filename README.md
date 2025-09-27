1. The crater code is a sample multifluid code written in F77 for a two-phase (air-water) model of a crater cavity simulation in 2D. The model system is the 5-equation transport equation with geometrical and gravitational source terms. The numerical method is the semi-discrete wave propagation method with THINC interface sharpening reconstruction.

2. The Nonlinearity_BBMH code is a code that accompanies the paper:  Hyperbolic approximation of the BBM equation, Nonlinearity, 33, 5477-5509, 2020.

3. The crater_clawpack is a mapped-grid version of the clawpack with a sample application to a model cavity problem. The model equation is the five-equation transport
   model. The equation of state is the stiffened gas equation of state

4. The carter_clawpack1 is a uniform Cartesian grid code with the same application as in the mapped grid code in 3.
