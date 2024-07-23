import os
import numpy as np

rundir = os.getcwd()
start = rundir+"/start.in"
run = rundir+"/run.in"
cparam = "./src/cparam.local"


def s_param():
    f = open(start, 'r')
    param = {}

    for line in f:
        line = line.strip()
        if not line.startswith(("!", "#", "&", "/")):
            line_comment_separated = line.split('!')[0].strip()
            # print(line_comment_separated.rstrip(',').rstrip('.'))

            if line_comment_separated.startswith("xyz_units"):
                units = line_comment_separated.split(' = ')[1].strip()
                conv_f = units.split(',')
                conv_fx = conv_f[0].strip()
                conv_fy = conv_f[1].strip()
                conv_fz = conv_f[2].strip()
                conv_f = [conv_fx, conv_fx, conv_fx]
                for ii in range(len(conv_f)):
                    if conv_f[ii]=="'pi'":
                        conv_f[ii] = np.pi
                    else:
                        ii = 1.0
                # print(conv_f[1])
            if line_comment_separated.startswith("xyz0"):
                xyz0 = line_comment_separated.split(' = ')[1].strip()
                # print(xyz0)
                xyz0 = xyz0.split(',')
                # print(xyz0)
                lx0 = float(xyz0[0].strip())*conv_f[0]
                ly0 = float(xyz0[1].strip())*conv_f[1]
                lz0 = float(xyz0[2].strip())*conv_f[2]
                # print(lx0, ly0, lz0)
                param.update({'lx0': lx0})
                param.update({'ly0': ly0})
                param.update({'lz0': lz0})

            elif line_comment_separated.startswith("Lxyz"):
                Lxyz = line_comment_separated.split(' = ')[1].strip()
                # print(Lxyz)
                Lxyz = Lxyz.split(',')
                # print(Lxyz)
                lx = float(Lxyz[0].strip())*conv_f[0]
                ly = float(Lxyz[1].strip())*conv_f[1]
                lz = float(Lxyz[2].strip())*conv_f[2]
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


def r_param():
    f = open(run, 'r')
    param = {}

    for line in f:
        line = line.strip()
        if not line.startswith(("!", "#", "&", "/")):
            line_comment_separated = line.split('!')[0].strip()
            # print(line_comment_separated.rstrip(',').rstrip('.'))
            if len(line_comment_separated)==0:
                continue
            else:
                line_comment_separated = line_comment_separated.rstrip(',').rstrip('.')
                line_comma_separated = line_comment_separated.split(',')
                for k_v in range(len(line_comma_separated)):
                    key, value = line_comma_separated[k_v].split('=')
                    key = key.strip()
                    value = value.strip()
                    # print(key)
                    # print(value)
                    try:
                        value = float(value)
                        param.update({key: value})
                    except ValueError:
                        if value[0]=="'":
                            value = value[1:-1]
                        param.update({key: value})
    return param


def c_param():
    f = open(cparam, 'r')
    param = {}

    for line in f:
        line = line.strip()
        if not line.startswith(("!", "#", "&", "/")):
            line_comment_separated = line.split('!')[0].strip()
            # print(line_comment_separated.rstrip(',').rstrip('.'))
            if len(line_comment_separated)==0:
                continue
            else:
                line_ratio_separated = line_comment_separated.split(' :: ')[1].strip()
                # print(line_ratio_separated.rstrip(',').rstrip('.'))

                line_ratio_separated = line_ratio_separated.rstrip(',').rstrip('.')
                line_comma_separated = line_ratio_separated.split(',')
                for k_v in range(len(line_comma_separated)):
                    key, value = line_comma_separated[k_v].split('=')
                    key = key.strip()
                    value = value.strip()
                    try:
                        value = float(value)
                        param.update({key: value})
                    except ValueError:
                        param.update({key: value})

    if param['nygrid']=='nxgrid':
        param['nygrid'] = param['nxgrid']

    if param['nzgrid']=='nxgrid':
        param['nzgrid'] = param['nxgrid']

    if param['nzgrid']=='nygrid':
        param['nzgrid'] = param['nygrid']

    return param


sparam = s_param()
rparam = r_param()
cparam = c_param()


nz = int(cparam['nzgrid'])    #nzgrid
depth = sparam['lz']   #Lz
# if sparam['']
z2 = sparam['lz0']+depth    #hight in z direction above Lref=0
z1 = sparam['lz0'] #depth in z direction below Lref=0

try:
    g = -rparam['gravz'] #acceleration due to gravity
except KeyError:
    try:
        g = -sparam['gravz']
    except:
        g = 1


try:
    gamma = rparam['gamma']
except:
    try:
        gamma = sparam['gamma']
    except:
        gamma = 5/3

print('z/L[0]:', round(z1,5), '\nz/L[-1]:', round(z2,5))
print('height of the box:', round(depth,5))
print('gravz:', -g)

w = rparam['wcool']    #parameter to control the sharpness of discontinuity

if str(rparam['cooltype'])=='rho_cs2':
    print("cooltype:", rparam['cooltype'])
    cs2max = rparam['cs2cool']    #cs^2 in the upper layer
    print("cs2max", cs2max)
    cs2base = 1.0
    print("cs2min", cs2base)
elif rparam['cooltype']=='two-layer':
    print("cooltype:", rparam['cooltype'])
    cs2max = rparam['cs2cool']    #cs^2 in the upper layer
    print("cs2max", cs2max)
    cs2base = rparam['cs2cool2']    # NOTE:
    print("cs2min", cs2base)


z = np.linspace(z1,z2,nz)   #grid
zz = np.zeros(np.size(z))   #array to save the integral values

cs2 = cs2base + 0.5*(1+np.tanh(z/w))*(cs2max-cs2base)

f = g/cs2
dz = z[1] - z[0]

for i in range(1, len(z)):
    # print(i)
    zz[i] = zz[i-1]+0.5*(f[i]+f[i-1])
zz = zz*dz  #value

lncs2 = np.log(cs2) #value to scale # NOTE: is it choosen arbitarily
lnrho = -gamma*(zz-depth)-lncs2 #log of \rho
ss = (gamma-1)*(zz-depth)+lncs2 #entropy

if 'cp' in sparam:
    cp = sparam['cp']
elif 'cp' in rparam:
    cp = rparam['cp']
else:
    cp = 1.0
# print('c_p:', cp)
R = (1-1/gamma)*cp
# print('R:', R)
TT = cs2/(gamma*R)

arr = np.empty([nz,3])
arr[:,0] = z
arr[:,1] = lnrho
arr[:,2] = ss

fmt = '%1.6f', '%1.6f', '%1.6f'

np.savetxt('stratification.dat', arr, delimiter=' ', fmt=fmt)

print('stratification.dat file has been created.')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, axs = plt.subplots(2,1, sharex=True, figsize=(5,3))
axs[0].plot(z, lnrho)
# axs[1].plot(z, ss)
axs[1].plot(z, TT)
axs[1].set_xlim(z[0], z[-1])
axs[1].set_xlabel('x')
axs[0].set_ylabel('lnrho')
axs[1].set_ylabel('TT')
plt.tight_layout()
plt.savefig('plots/initial_strati.png')