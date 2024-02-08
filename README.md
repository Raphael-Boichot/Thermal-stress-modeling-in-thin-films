# Thermal stress modeling in thin films

This MATLAB set of codes is intended to predict the stress evolution in a system comprising a substrate and one layer or several layers of coatings and subjected to a variation of temperature (typically from fabrication high temperature to room temperature). The Optimization function to solve this problem is particularly tricky to minimize so the reader would not be surprised that the minimization method is non canonic (it couples a Nelder-Meads algorithm to a Genetic Algorithm). It easily falls into local minima so the code is reloaded automatically in case a local minima is detected (typically a minimum not close enough to zero). 

In this model, the "film" can be as thick as the "substrate", there is no limit for their ratio of thickness. The film is just supposed to be on top of the substrate for ease of understanding.

## Single thin film on susbtrate ("I'm Too Young to Die")

Strain calculation in the multilayer follows the classical beam theory and the formalism of [Hsueh et al.](https://doi.org/10.1016/S0040-6090(02)00699-5). It includes the effect of island growth following the model of [Wu et al.](https://doi.org/10.1557/PROC-0892-FF26-01) for the particular case of AlN growth on sapphire. It can be adapted to any other couple of materials, non ceramic included. The code returns the position of neutral axis, the final stress and the radius of curvature of the system at mechanical equilibrium.

How to use it : 
- **[Minimize_GA.m](https://github.com/Raphael-Boichot/Thermal-stress-modeling-in-thin-films/blob/main/Codes%20single%20layer/Minimize_GA.m)** : main file where geometric parameters of subtrate and thin film, growth strain in the coating and temperature can be modified. The thermal expansion coefficients for both materials are entered as polynomials for maximum accuracy. All mechanical setting parameter and properties are [here](https://github.com/Raphael-Boichot/Thermal-stress-modeling-in-thin-films/blob/f6da0caecd306a09ab7a5f4855436a23501f7c9d/Codes%20single%20layer/Minimize_GA.m#L20).
- **[Contrainte_multicouches.m](https://github.com/Raphael-Boichot/Thermal-stress-modeling-in-thin-films/blob/main/Codes%20single%20layer/contrainte_multicouches.m)**: mechanical modeling itself based on Hsueh et al. The model consists in zeroing independently 3 beefy integrals at the same time (namely the force due to the uniform strain, force due to the bending strain and the sum of bending moments) by varying the total stress, the redius of curvature and the position of neutral axis. There is an obvious local minima where the curvature is inversed (the solution is a bell shaped surface) so the code systematically searches for better solutions at the opposite of a converged one.
- **[GA.m](https://github.com/Raphael-Boichot/Thermal-stress-modeling-in-thin-films/blob/main/Codes%20single%20layer/GA.m) and [croisement_GA.m](https://github.com/Raphael-Boichot/Thermal-stress-modeling-in-thin-films/blob/main/Codes%20single%20layer/croisement_GA.m)**: calculation of Genetic algorithm coupled to a Nelder-Mead algorithm. Typically each candidate population is first subjected to a gradient free optimization (Nelder-Mead) then to a global optimization (Genetic Algorithm). Normally the inverse (global then local optimization) is what to do but not in this case. If you do not understand what I say, do not modify any parameter into these functions.

Result should be close to the simplified Stoney equation in case the film is infinitely thin compared to the substrate but gives curvature and stress in any case. There is no approximation in this method. To my knowledge, the codes always converges to the global optimum solution as it is written.

## Multi-layer thin film on susbtrate subjected to thermal cycling and creeping ("Hurt Me Plenty")

To cite this work :
[Danying Chen, Alexandre Crisci, Raphaël Boichot, Johann Colas, Ludovic Charpentier, Marianne Balat-Pichelin, Michel Pons, Frédéric Mercier,
Modeling multilayer coating systems in solar receivers, Surface and Coatings Technology, Volume 399, 2020, 126102, ISSN 0257-8972, https://doi.org/10.1016/j.surfcoat.2020.126102.](https://www.sciencedirect.com/science/article/abs/pii/S0257897220307714)

Details are given in the [preprint](https://github.com/Raphael-Boichot/Thermal-stress-modeling-in-thin-films/blob/main/Codes%20multilayers%20with%20creeping/Documentation/articleCSP.pdf) and the [supplementary materials](https://github.com/Raphael-Boichot/Thermal-stress-modeling-in-thin-films/blob/main/Codes%20multilayers%20with%20creeping/Documentation/articleCSP_Supplementary%20Materials.pdf). Based on the same formalism and the same method of resolution as **Single thin film on substrate**.

How to use it :
- **[Main_thermal_cycling.m](https://github.com/Raphael-Boichot/Thermal-stress-modeling-in-thin-films/blob/main/Codes%20multilayers%20with%20creeping/Main_thermal_cycling.m)** : main file where geometric parameters of subtrate and thin film, thermal boundary conditions (T and heat flux), growth strain in the coating and temperature, cycle (holding time and cooling temperature can be modified. 
- **[cycle_GA.m](https://github.com/Raphael-Boichot/Thermal-stress-modeling-in-thin-films/blob/main/Codes%20multilayers%20with%20creeping/cycle_GA.m)** contains the definition of the thermal expansion coefficents (entered as polynomials for maximum accuracy). cycle_GA is used for calculating the stress evolution during thermal cycling.
- **[Contrainte_multicouche_Hsueh_4.m](https://github.com/Raphael-Boichot/Thermal-stress-modeling-in-thin-films/blob/main/Codes%20multilayers%20with%20creeping/Contrainte_multicouche_Hsueh_4couches.m)** contains the mechanical modeling itself based on Hsueh et al. The model consists in zeroing independently 3 beefy integrals at the same time (namely the force due to the uniform strain, force due to the bending strain and the sum of bending moments) by varying the total stress, the redius of curvature and the position of neutral axis. There is an obvious local minima where the curvature is inversed (the solution is a bell shaped surface) so the code systematically searches for better solutions at the opposite of a converged one.
- **[minimize_GA_TZM_1.m](https://github.com/Raphael-Boichot/Thermal-stress-modeling-in-thin-films/blob/main/Codes%20multilayers%20with%20creeping/minimize_GA_TZM_1.m)** is used for calculating initial stress distribution in the system (just after deposition)
- **[GA.m](https://github.com/Raphael-Boichot/Thermal-stress-modeling-in-thin-films/blob/main/Codes%20multilayers%20with%20creeping/GA.m) and [croisement_GA.m](https://github.com/Raphael-Boichot/Thermal-stress-modeling-in-thin-films/blob/main/Codes%20multilayers%20with%20creeping/croisement_GA.m)**: calculation of Genetic algorithm coupled to a Nelder-Mead algorithm. Typically each candidate population is first subjected to a gradient free optimization (Nelder-Mead) then to a global optimization (Genetic Algorithm). If you do not understand what I say, do not modify any parameter into these functions.

## Multi-layer thin film on susbtrate subjected to thermal cycling, creeping and oxydation ("Ultra violence")

% 1-code 'n_cycle' & 'cycle_GA' & 'Contrainte_multicouche_Hsueh_4' couche are used for claculating the stress evolution during thermal cycling
% 
% 2-code 'minimize_GA_TZM_1' is used for claculating initial stress distribution in the system (just after deposition)
% 
% 3-files 'AlN_TZM_AlN_q10000_xx_CTEofT' are obtained using instantaneous CTE for claculations
%   also, we used a deltat of 1/12 (h) and a node number of 5000 for the substrate during these claculations
% 
% 4-files 'AlN_TZM_AlN_q10000' are obtained using constant CTE
%   also, we used a deltat of 1/12 (h) and a node number of 5000 for the substrate during these claculations
% 
% 5-files 'AlN_TZM_AlN_maille50000_q10000_Tambxx_xx': use 25 °C or 20 °C as ambiant temperature in the code (for the others without Tamb mentioned in the file name, we use 25 °C as ambiant temperature)
%   also, we used a deltat of 1/6 (h) and a node number of 50000 for the substrate during these claculations
% 
% 5-files 'AlN_TZM_AlN_maille50000_q10000_xx_CTEofT': use instantaneous CTE for claculations
%   also, we used a deltat of 1/6 (h) and a node number of 50000 for the substrate during these claculations
