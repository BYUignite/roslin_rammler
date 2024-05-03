# Roslin Rammler distribution fit

* Find the Roslin Rammler distribution for an experimental cumulative size distribution.
    * Given the particle sizes (diameters) and the cumulative distribution fractions.
    * Fit the data to find the R-R parameters
* Given the R-R distribution, select a desired number N of particle sizes to consider.
    * Compute 2N moments of the distribution and invert these moments to find corresponding quadrature weights and abscissas.
        * The abscissas are the particle diamters.
        * The weights are the number of particles.
    * Scale the weights by a constant factor so that the total particle mass is some desired value: mcoal.

[Link to notebook file](https://nbviewer.org/github/BYUignite/roslin_rammler/blob/master/Roslin_Rammler.ipynb)
