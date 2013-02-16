import matplotlib.pyplot as plt
import numpy as np

amb = np.loadtxt("x")
ele = np.loadtxt("y")

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

ax1.plot(amb[:,0], "b.", label="Amb", alpha=0.4)
ax1.plot(ele[:,0], "r.", label="Ele", alpha=0.4)
ax1.legend(loc="best", ncol=1)
ax1.set_ylabel("GPP (gC m$^{-2}$ d$^{-1}$)")


ax2.plot(ele[:,0]/amb[:,0], "g-")
ax2.set_ylabel("Response (Ele/Amb)")
ax2.set_xlabel("Doy since 1996")
plt.show()