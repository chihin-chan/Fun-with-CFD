# Fun-with-CFD
I created these codes to familiarise myself with the nature of hyperbolic PDEs. Most of the numerical schemes can be found in the textbook  'Numerical methods for conservation laws' by Randall Leveque

In the code, the user is allowed to choose various numerical schemes and initial conditions.

Numerical Schemes
1. Lax-Friedrichs
2. Lax-Wendroff
3. Ritchymyer
4. MarCormack

Initial Condition
1. Gaussian
2. N Function
3. Inverted N Function
4. Staircase

From the code, I observed that the CFL condition of < 1 is not sufficient as the Burger's equation is nonlinear. An intial condition of N function with Lax-Wendroff numerical scheme will also yield a non-physical result which can be fixed by the Ritchmyer scheme.
