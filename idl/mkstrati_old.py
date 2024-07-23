import numpy as np

nz = 256    #nzgrid
depth = np.pi   #Lz
z2 = np.pi/10    #hight in z direction above Lref=0
z1 = z2 - depth #depth in z direction below Lref=0
g = 1.0 #acceleration due to gravity
gamma = 5/3
w = 0.06    #parameter to control the sharpness of discontinuity
cs2max = 4.0    #cs^2 in the upper layer
csbase = 1.0    # NOTE:

z = np.linspace(z1,z2,nz)   #grid
zz = np.zeros(np.size(z))   #array to save the integral values

cs2 = csbase + 0.5*(1+np.tanh(z/w))*(cs2max-1)

f = g/cs2
dz = z[1] - z[0]

for i in range(1, len(z)):
    # print(i)
    zz[i] = zz[i-1]+0.5*(f[i]+f[i-1])
zz = zz*dz  #value

lncs2 = np.log(cs2) #value to scale # NOTE: is it choosen arbitarily
lnrho = -gamma*(zz-depth)-lncs2 #log of \rho
ss = (gamma-1)*zz+lncs2 #entropy

arr = np.empty([nz,3])
arr[:,0] = z
arr[:,1] = lnrho
arr[:,2] = ss

np.savetxt('stratification.dat', arr, delimiter=' ')