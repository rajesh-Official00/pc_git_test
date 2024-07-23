import os
import numpy as np

rundir = os.getcwd()
start = rundir+"/start.in"


def s_param():
    f = open(start, 'r')
    param = {}

    for line in f:
        line = line.strip()
        if not line.startswith(("!", "#", "&", "/")):
            line_comment_separated = line.split('!')[0].strip()
            # print(line_comment_separated.rstrip(',').rstrip('.'))

            if line_comment_separated.startswith("xyz0"):
                xyz0 = line_comment_separated.split(' = ')[1].strip()
                # print(xyz0)
                xyz0 = xyz0.split(',')
                # print(xyz0)
                lx0 = float(xyz0[0].strip())
                ly0 = float(xyz0[1].strip())
                lz0 = float(xyz0[2].strip())
                # print(lx0, ly0, lz0)
                param.update({'lx0': lx0})
                param.update({'ly0': ly0})
                param.update({'lz0': lz0})

            elif line_comment_separated.startswith("Lxyz"):
                Lxyz = line_comment_separated.split(' = ')[1].strip()
                # print(Lxyz)
                Lxyz = Lxyz.split(',')
                # print(Lxyz)
                lx = float(Lxyz[0].strip())
                ly = float(Lxyz[1].strip())
                lz = float(Lxyz[2].strip())
                # print(lx, ly, lz)
                param.update({'lx': lx})
                param.update({'ly': ly})
                param.update({'lz': lz})

            elif line_comment_separated.startswith(("gamma", "gravz")):
                line_comment_separated = line_comment_separated.rstrip(',').rstrip('.')
                line_comma_separated = line_comment_separated.split(',')
                for k_v in range(len(line_comma_separated)):
                    key, value = line_comma_separated[k_v].split('=')
                    key = key.strip()
                    value = value.strip()
                    try:
                        value = float(value)
                        param.update({key: value})
                    except ValueError:
                        param.update({key: value})
    return param



def k_vec(dkx, dky, dkz, k1, k2, kmax):
    # dkx=1.; dky=1.; dkz=1.
    ex=1.; ey=1.; ez=1.
    # k1=4.2; k2=5.5
    # kmax=10.

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



sparam = s_param()

lx = sparam['lx']
ly = sparam['ly']
lz = sparam['lz']

dkx = np.round(2*np.pi/lx, 6)
dky = np.round(2*np.pi/ly, 6)
dkz = np.round(2*np.pi/lz, 6)

print(dkx, dky, dkz)

k_vec(dkx, dky, dkz, 7.5, 8.5, 14)