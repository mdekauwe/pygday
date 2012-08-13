import numpy as np
import matplotlib.pyplot as plt

amb = np.loadtxt("x")
ele = np.loadtxt("xx")
effect = []
for yr in np.arange(1998, 2008):
    a = np.sum(amb[np.where(amb[:,0]==yr)][:,1])
    e = np.sum(ele[np.where(ele[:,0]==yr)][:,1])
    effect.append(e/a)
print np.mean(effect)
plt.plot(effect)
plt.show()