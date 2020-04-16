# -*- coding: utf-8 -*-
"""
NAME: m2hPa(z)

DESCRIPTION: estimate atmospheric pressure at given altitude above mean sea level
             Air pressure at sea level is considered 1013.25 hPa
             www.weather.gov/epz/wxcalc_pressurealtitude

INPUT:
    z, meters above sea level (m)
    
OUTPUT:
    float, atmospheric pressure (hPa)

@author: Daniele Zannoni
Created on Sun Mar 29 21:50:40 2020
"""
def  m2hPa(m):
    return (1-(m/44307.69396))**(1/0.190284)*1013.25


"""
NAME: hPa2m(hPa)

DESCRIPTION: estimate altitude at given atmospheric pressure above mean sea level
             Air pressure at sea level is considered 1013.25 hPa
             www.weather.gov/epz/wxcalc_pressurealtitude

INPUT:
    hPa, atmospheric pressure (hPa)
    
OUTPUT:
    float, meters above sea level (m)

@author: Daniele Zannoni
Created on Sun Mar 29 21:50:40 2020
"""
def hPa2m(hPa):
    return (1-(hPa/1013.25)**0.190284)*44307.69396

"""
NAME: rhoAir(z, T)

DESCRIPTION: estimate air density at given altitude and temperature above 
             mean sea level. Air pressure at sea level is considered 1013.25 hPa

INPUT:
    z, meters above sea level (m)
    T, air temperature at altitude z (°C)
    
OUTPUT:
    float, air density (kg m^-3)

@author: Daniele Zannoni
Created on Sun Mar 29 21:50:40 2020
"""
def rhoAir(z, T):
    # Universal gas constant for dry air (J kg^-1 K^-1)
    # https://www.skybrary.aero/index.php/Density_Altitude
    Rair = 287.058
    # Convert °C to K
    T = T + 273.15
    # Estimate air pressure at given height
    airP = (1-(z/44307.69396))**(1/0.190284)*1013.25
    # Return air density 
    return 100*airP/(Rair*T)

"""
NAME: SatVapPress(T)

DESCRIPTION: estimate saturation water vapor pressure (Pa) at given temperature (°C)
using the Buck formula which has the lowest error in temperature range 0-35°C
References: Buck, A. L. (1981), "New equations for computing vapor pressure 
and enhancement factor", J. Appl. Meteorol., 20: 1527–1532

INPUT:
    T, temperature (°C)
    
OUTPUT:
    float, water vapor pressure (hPa)

@author: Daniele Zannoni
Created on Sun Mar 29 21:50:40 2020
"""

def SatVapPress(T):
    from math import exp # required exp
    return 6.1121*exp((18.678-T/234.5)*(T/(257.14+T)))