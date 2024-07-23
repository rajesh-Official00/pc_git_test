import numpy as np
import matplotlib.pyplot as plt
import pencil as pc

#reading the xy_aver file
xya=pc.read.aver(plane_list=['xy'])
bxmz=xya.xy.bxmz
bymz=xya.xy.bymz
by = [bxmz, bymz]

#reading the grid
grid = pc.read.grid(trim=True)
z = grid.z

fig, axs = plt.subplots(2,1, figsize=(10,8), sharex=True)
plt.rcParams.update({'font.size': 16})
plt.rcParams['text.usetex'] = True
plt.rc("figure.subplot", left=0.2)
plt.rc("figure.subplot", right=0.95)
plt.rc("figure.subplot", bottom=0.15)
plt.rc("figure.subplot", top=0.90)

for i in range(len(axs.flat)):
    ims = axs[i].contour(np.transpose(by[i]),50,cmap='bwr')

#axs[0].contour(np.transpose(bxmz),50,cmap='bwr')
axs[0].set_xlabel(r'$t$')
axs[0].set_ylabel(r'$z$')
axs[0].set_title(r'$\langle B_x \rangle _{xy}$')
axs[0].set_xlim(0,)

#axs[1].contour(np.transpose(bymz),50,cmap='bwr')
axs[1].set_xlabel(r'$t$')
axs[1].set_ylabel(r'$z$')
axs[1].set_title(r'$\langle B_y \rangle _{xy}$')
axs[1].set_xlim(0,)

#cbar_ax = fig.add_axes([1.05, 0.15, 0.05, 0.7])
#fig.colorbar(ims, cax=cbar_ax)
plt.tight_layout()
#plt.savefig('ave.jpg')

plt.show()
