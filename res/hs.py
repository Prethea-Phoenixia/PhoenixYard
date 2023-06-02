"""
Hirschfelder-Sherman calculation
additive approximation for propellant thermochemical properties,
by tallying propellant heat of explosion and calculating the
deviaiton to the reference temperature of 2500K, assuming 
constant heat capacity between 2500K and 3000K, propellant isochoric
adiabatic flame temperature is then calculated. 

Combustion product species are assumed to involve CO, H2O, CO2, H2,
N2, HCl , and linear equations are solved to estimate product gas 
volume (the inverse of molar mass), and adiabatic ratios.

Results are within few percent of accurate thermalchemical calculation
in the range of 2000K to 4000K. Should not be applied to cases >4000K
where dissociation of of free radicals like -H, -OH and -Cl would happen.
Should not be applied to propellants generating substantial amount of 
consdensed exhaust.


Significant amount of propellant research nowadays are
published with incomplete data to be used in the inetrnal
ballistic formulation that is present in this program. Ideally
accurate thermochemcial programs like BLAKE should be consulted.
In the absence of more sophisticated programs, the Hirschfelder
-Sherman is an approximate, and simple method appropriate for
propellant screening and program guidance purposes.
"""
