import numpy as np
import matplotlib.pyplot as plt
import pencil as pc

ts= pc.read.ts()

plt.plot(ts.t, ts.brms)
plt.yscale('log')
#plt.ylim(0,100)
plt.show()
