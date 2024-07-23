#!/usr/bin/python
# script for plotting different diagrams (timeseries, k-omega, butterfly)
# diagram using module different modules which in turn use module modes

import os
import sys
import pencil as pc
import matplotlib.pyplot as plt

# location of modules
modulePath = '/scratch/cmondal/pc_data/dynamo_modes/test_dir/my_scripts/'
sys.path.append(modulePath)

# importing my modules
import statevariables
import rms
import butterfly
import kom

plt.rcParams.update({'font.size': 14})
plt.rcParams['text.usetex'] = True

# location of diagrams
path = './test/'

pathexists = os.path.exists(path)
if pathexists:
    print(path + ' exists.')
else:
    print(path + ' does not exist.')
    os.makedirs(path)
    print(path + ' is created.')

# reading the simulation object and other files
sim = pc.get_sim(quiet=True)
ts = pc.read.ts()
yaver = pc.read.aver(plane_list='y')
xyaver = pc.read.aver(plane_list=['xy'])

# # creating an StateVariables object
# img_StateVariables = statevariables.StateVariables(path, sim, ts, xyaver, yaver, ini=False, dyn=True)
# # plotting butterfly diagram
# img_StateVariables.plot(8,5)


# # creating an Rms object
# img_Rms = rms.Rms(path, sim, ts, xyaver, yaver, ini=False, dyn=True)
# # plotting butterfly diagram
# img_Rms.plot(8,5)


# # creating an Butterfly object
# img_Bfly = butterfly.Butterfly(path, sim, ts, xyaver, yaver, ini=False, dyn=True)
# # plotting butterfly diagram
# img_Bfly.plot(ticks=[-.3, -.2, -.1, 0, .1, .2, .3])


# creating an KOm object
img_KOm = kom.KOm(path, 150, 670, 230, sim, ts, xyaver, yaver, ini=False, dyn=True)
# plotting k-omega diagram
img_KOm.plot(cmap='hot_r')