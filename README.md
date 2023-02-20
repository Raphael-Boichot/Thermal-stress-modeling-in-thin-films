# Thermal stress modeling in thin films

This MATLAB set of codes is intended to predict the stress evolution in a system comprising a substrate and one layer or several layers of coating and subject to a variation of temperature (typically from fabrication temperature to room temperature).  The Optimization function to solve this problem is particularly tricky to minimize so the reader would not be surprised that the minimization method is non canonic (it couples a Nelder-Meads alogorithm to a Genetic Algorithm). It easily falls into local minima so the code can be reloaded automatically in case a local minima is detected (typically a minimum not close enough to zero).

**Single thin film on susbtrate**

Strain calculation in the multilayer follows the classical beam theory and the formalism of [Hsueh et al.](https://doi.org/10.1016/S0040-6090(02)00699-5). It includes the effect of island growth following the model of [Wu et al.](https://doi.org/10.1557/PROC-0892-FF26-01) for the particular case of AlN growth on sapphire. It can be adapted to any other couple of materials, non ceramic included. The code returns the position of neutral axis, the final stress and the radius of curvature of the system at mechanical equilibrium.

How to use it : 
- **Minimize_GA** : main file where geometric parameters of subtrate and thin film, growth strain in the coating and temperature can be modified.
- **Contrainte_multicouche_BOICHOT**: mechanical modeling itself based on Hsueh et al.
- **GA and croisement GA**: calculation of Genetic algorithm coupled to a Neledr-Mead algorithm.

**Multi-layer thin film on susbtrate

TO DO
