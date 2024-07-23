import numpy as np

dkx=1.; dky=1.; dkz=1.
ex=1.; ey=1.; ez=1.
k1=3.6; k2=4.5
kmax=10.

kav = 0
# kmaxz = kmax

if kmax<k2:
    print('Warning: non-spherical region in k-space')

kx_range = np.arange(-kmax, kmax, dkx)
ky_range = np.arange(-kmax, kmax, dky)
kz_range = np.arange(-kmax, kmax, dkz)

kkx = []
kky = []
kkz = []

# i = 0
for kx in kx_range:
    for ky in ky_range:
        for kz in kz_range:
            k = np.sqrt(kx**2+ky**2+kz**2)
            kref = np.sqrt((kx/ex)**2+(ky/ey)**2+(kz/ez)**2)
            if kref>k1 and kref<k2:
                kav = kav+k
                kkx.append(kx)
                kky.append(ky)
                kkz.append(kz)

n = len(kkx)
kav = round(kav/n,5)

print('Writing', n, 'wave vectors; kav = ', kav)

kkx = np.array(kkx)
kky = np.array(kky)
kkz = np.array(kkz)

f=open('k.dat','w')

x = 6
# y = 6

f.write(str(n)+"      "+ str(kav)+ "\n")
f.write("\n")

for i in range(0,len(kkx),x):
    np.savetxt(f, kkx[i:i+x], fmt='%7.4f', delimiter="    ", newline=", ")
    f.write("\n")
f.write("\n")

for i in range(0,len(kky),x):
    np.savetxt(f, kky[i:i+x], fmt='%7.4f', newline=", ")
    f.write("\n")
f.write("\n")

for i in range(0,len(kkz),x):
    np.savetxt(f, kkz[i:i+x], fmt='%7.4f', newline=", ")
    f.write("\n")

f.close()

print('check for isotropy: <k>=%8.5f, %8.5f, %8.5f' %(np.average(kkx), np.average(kky), np.average(kkz)))
print('check for isotropy: <k^2>= %8.5f, %8.5f, %8.5f' %(np.average(kkx**2), np.average(kky**2), np.average(kkz**2)))
