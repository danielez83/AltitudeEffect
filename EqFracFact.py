# -*- coding: utf-8 -*-
"""
NAMES: alpha18(T), alpha2(T)

DESCRIPTION: calculate equilibrium fractionation factor for Oxigen18/Oxygen16 and
Deuterium/Hydrogen ratios on liquid/vapor system at given temperature (°C) 
between -5°C and +100°C with the formula of Majoube (1971)
below -5°C with the formula of Ellehoj et al. (2013)
REFERENCES: 
    Majoube, Michel. "Fractionnement en oxygene 18 et en deuterium entre l’eau et sa vapeur." 
Journal de Chimie Physique 68 (1971): 1423-1436.
    Ellehoj, M. D., et al. "Ice‐vapor equilibrium fractionation factor of hydrogen and oxygen isotopes: 
Experimental investigations and implications for stable water isotope studies." 
Rapid Communications in Mass Spectrometry 27.19 (2013): 2149-2158.


INPUT:
    T, temperature (°C)
    
OUTPUT:
    float, liquid/vapor fractionation factor (-)

@author: Daniele Zannoni
Created on Sun Mar 29 21:50:40 2020
"""

def alpha18(T):
    from math import exp # required exp
    T = T+273.15 # Convert °C to K
    if T>=268.16:
        return exp((1137/(T**2))-(0.4156/T)-0.0020667)
    else:
        return exp((8312.58/(T**2))-(49.192/T)+0.0831)

def alpha2(T):
    from math import exp # required exp
    T = T+273.15 # Convert °C to K
    if T>=268.16:    
        return exp((24844/(T**2))-(76.248/T)+(0.052612))
    else:
        return exp((48888/(T**2))-(203.10/T)+0.2133)
