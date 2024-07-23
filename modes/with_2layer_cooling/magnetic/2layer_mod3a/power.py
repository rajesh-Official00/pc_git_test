import numpy as np
import matplotlib.pyplot as plt
import pencil as pc

p = pc.read.power()

uz = np.fft.fft(p.uz_xy, axis=0)
plt.imshow(np.abs(uz[:1000,0,:]))
plt.colorbar()
plt.ylim(0,100)
plt.show()
