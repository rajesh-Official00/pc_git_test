import pencil as pc
import matplotlib.pylab as plt
import matplotlib.colors as colors
import numpy as np

xya=pc.read.aver(plane_list=['xy'])
bxmz=xya.xy.bxmz
bymz=xya.xy.bymz

plt.contour(np.transpose(bxmz),50,cmap='bwr')
plt.xlabel(r'$t$')
plt.ylabel(r'$z$')

plt.tight_layout()

plt.show()
