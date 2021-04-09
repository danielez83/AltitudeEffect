# -*- coding: utf-8 -*-
"""
Description:    Estimate isotopic lapse rate (altitude effect) by using a simple 
                Rayleigh distillation model
Author:         Daniele Zannoni
Date:           08/04/2020

Requirements:   EqFracFact.py (equilibrium fractionation factors estimation)
                MetFun.py (set of functions to calculate/convert meteorological
                parameters)

Notes:          Mean Lapse rate from Minder et al. (2010)
                Minder, Justin R., Philip W. Mote, and Jessica D. Lundquist. 
                "Surface temperature lapse rates over complex terrain: Lessons 
                from the Cascade Mountains." Journal of Geophysical Research: 
                Atmospheres 115.D14 (2010).
                
                VSMOW ratios from wikipedia
                https://en.wikipedia.org/wiki/Vienna_Standard_Mean_Ocean_Water
"""

import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import numpy as np
import math

# Customm functions
import EqFracFact
import MetFun

#%% Settings
# Simulation parameters and boundary conditions
MLR         = 5                 # °C/km, Moist adiabatic lapse rate 
MinZ        = 0                 # m, lowest altitude of the simulation
MaxZ        = 4000              # m, top level of the simulation
Zres        = 10                # m, vertical resolution of the simulation

# Initial conditions of the simulation
Tz0         = 20                # °C, air temperature at z=MinZ m ASL
P0          = 1013.25           # hPa, atmospheric pressure at z=0 
RH0         = 85/100            # %, Water vapor concentration at initial conditions (mixing ratio) ~ 19000 # ppmv
d18O0       = -12               # ‰, Water vapor initial d18O 
dD0         = 10 + (8*d18O0)    # ‰, Water vapor initial dD (d-excess = 10‰)

# Other parameters
VSMOW1816   = 2005.20*1e-6      # O-18/O-16 absolute ratio in VSMOW
VSMOW21     = 155.76*1e-6       # H-2/H-1 absolute ratio in VSMOW

# Vectors and matrixes memory allocation
Zvector                 = np.linspace(MinZ, MaxZ, np.int((MaxZ - MinZ)/Zres))   # Altitude levels 
Tvector                 = Tz0 - (Zvector*MLR/1000)                              # Air temperature 
VaporConcentration      = np.zeros(shape=(len(Zvector),1))                      # Water vapor concentration 
VaporComposition        = np.zeros(shape=(len(Zvector), 2))                     # Water vapor composition. Column 1 d18O, Column 2 dD
VaporFraction           = np.zeros(shape=(len(Zvector),1))                      # Residual fraction of water vapor
PrecComposition         = np.zeros(shape=(len(Zvector), 2))                     # Precipitation composition. Column 1 d18O, Column 2 dD
RHvector                = np.zeros(shape=(len(Zvector),1))                      # Relative Humidity

#%% Run Simulation
# Assign initial conditions 
VaporConcentration[0]   = 1e6*(RH0)*MetFun.SatVapPress(Tz0)/P0
VaporFraction[0]        = 1
RHvector[0]             = RH0
VaporComposition[0,0]   = d18O0
VaporComposition[0,1]   = dD0
R1816_0_VAP             = ((d18O0*1e-3)+1)*VSMOW1816
R21_0_VAP               = ((dD0*1e-3)+1)*VSMOW21

# Lift the air mass
print("Lifting air mass from %.0f m to %.0f m...\n"%(MinZ, MaxZ))
for z in range(1, len(Zvector)):
    # RH at current z level
    RHvector[z] = (((VaporConcentration[z-1]/1e6)*P0))/(MetFun.SatVapPress(Tvector[z]))
    # Condensation?
    if RHvector[z]>1.05: # just little higher than 100% to account for error in the sat press formula
        # New mixing ratio
        VaporConcentration[z] = 1e6*MetFun.SatVapPress(Tvector[z])/MetFun.m2hPa(z)
        # Estimate fraction of residual water vapor
        VaporFraction[z] = VaporConcentration[z]/VaporConcentration[0]
        # Estimate isotopic composition of water vapor
        R1816_VAP = R1816_0_VAP*VaporFraction[z]**(EqFracFact.alpha18(Tvector[z])-1)
        R21_VAP = R21_0_VAP*VaporFraction[z]**(EqFracFact.alpha18(Tvector[z])-1)
        # Convert to delta units
        VaporComposition[z, 0] = (R1816_VAP/VSMOW1816 - 1)*1e3
        VaporComposition[z, 1] = (R21_VAP/VSMOW21 - 1)*1e3
        # Estimate isotopic composition of precipitation by mass balance
        PrecComposition[z,0] = (VaporFraction[z-1]*VaporComposition[z-1, 0]-VaporFraction[z]*VaporComposition[z, 0])/(VaporFraction[z-1]-VaporFraction[z])
        PrecComposition[z,1] = (VaporFraction[z-1]*VaporComposition[z-1, 1]-VaporFraction[z]*VaporComposition[z, 1])/(VaporFraction[z-1]-VaporFraction[z])       
    else:
        VaporConcentration[z] = VaporConcentration[z-1]
        VaporFraction[z] = VaporFraction[z-1]
        VaporComposition[z, 0] = VaporComposition[z-1, 0]
        VaporComposition[z, 1] = VaporComposition[z-1, 1]
        PrecComposition[z, 0] = PrecComposition[z-1, 0]
        PrecComposition[z, 1] = PrecComposition[z-1, 1]

# Clean Precipitation array, i.e. remove all values (zeros) before the first condensation event
i = 0
while PrecComposition[i, 0] == 0:
    PrecComposition[i, :] = math.nan
    i = i + 1
print("End of simulation...")

#%% Sample randomly 100 precipitation and altitude values to build a linear model
rand_indexes = np.random.randint(i, len(Zvector)-1, 20)# np.random.random_integers(i, len(Zvector)-1, 20)
dummy_altitude = Zvector[rand_indexes]
dummy_prec_composition = PrecComposition[rand_indexes, 0]
# Calculate linear regression
regressor = LinearRegression()
regressor.fit(dummy_altitude.reshape(-1, 1), dummy_prec_composition.reshape(-1, 1))
predicted_composition = regressor.predict(dummy_altitude.reshape(-1, 1))
R2 = r2_score(dummy_prec_composition.reshape(-1, 1), predicted_composition)
print("Isotopic lapse rate estimated: %.2f ‰/100m" % (regressor.coef_*100))

#%% Make plots
# Prepare subplot figure
ax1=plt.subplot(2, 1, 1)
# Plot water vapor isotopic composition as a function of altitude
plt.plot(Zvector[:], VaporComposition[:, 0])
plt.plot(Zvector[:], PrecComposition[:, 0])
plt.ylabel(r'$\delta^{18}$O (‰)')
plt.legend(['Water Vapor', 'Precipitation'])
# Plot regression model
ax2=plt.subplot(2, 1, 2)
plt.scatter(dummy_altitude, dummy_prec_composition)
plt.plot(dummy_altitude, predicted_composition, color='red', linewidth=1, linestyle='dashed')
# Show model information
buff_txt = r"$\delta^{18}$O=%.4f*Altitude%.4f (R$^{2}$=%.3f)" % (regressor.coef_, regressor.intercept_, R2)
plt.text(0.5*max(dummy_altitude), max(dummy_prec_composition)-.5, buff_txt)
plt.xlabel('Altitude (m ASL)')
plt.ylabel(r'$\delta^{18}$O (‰)')
