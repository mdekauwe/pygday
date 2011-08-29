#!/usr/bin/env python

"""
Make a quick plot of the output file...

that's all folks.
"""
__author__  = "Martin De Kauwe"
__version__ = "1.0 (21.03.2011)"
__email__   = "mdekauwe@gmail.com"

import numpy as np
import sys
import matplotlib.pyplot as plt


mate = np.loadtxt("z")
#mate_water = np.loadtxt("zz")
bewdy = np.loadtxt("y")
assim_mod3 = np.loadtxt("x")

fig = plt.figure()
plt.rcParams.update({'legend.fontsize': 8})
plt.rcParams.update({'mathtext.default': 'regular'})
ax = fig.add_subplot(111)
ax.plot(assim_mod3, '.', label="Assimilation Model 3" )

ax.plot(bewdy, '.', label="BEWDY")
ax.plot(mate, '.', label="MATE")
#ax.plot(mate_water, '.', label="MATE Water")

ax.legend(numpoints=1, loc='best', shadow=True).draw_frame(True)
ax.set_xlabel("Day of year")
ax.set_ylabel("NPP (gCm2day)")
fig.savefig("/Users/mdekauwe/Desktop/npp_comp.png")

