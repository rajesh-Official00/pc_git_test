#!/usr/bin/python
# script for plotting different diagrams (timeseries, k-omega, butterfly)
# diagram using module different modules which in turn use module modes

import os
import sys
import pencil as pc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#!! location of modules
modulePath = '/scratch/cmondal/pc_data/dynamo_modes/test_dir/scripts/for_modes'
sys.path.append(modulePath)

#!! importing my modules
# import statevariables
# import rms
import butterfly
# import kom

plt.rcParams.update({'font.size': 14})
# plt.rcParams['text.usetex'] = True


#!! location of working directory
script_path = os.getcwd()
sim_path = os.path.abspath(os.path.join(script_path))#, os.pardir))
print('Simulation directory: ', sim_path)
data_dir = '/data'


print('importing data ...')
# reading the simulation object and other files
sim = pc.get_sim(path=sim_path, quiet=True)
ts = pc.read.ts(datadir=sim_path+data_dir, quiet=True)
yaver = pc.read.aver(datadir=sim_path+data_dir, simdir=sim_path, plane_list='y')
xyaver = pc.read.aver(datadir=sim_path+data_dir, simdir=sim_path, plane_list=['xy'])

# cursor up one line
sys.stdout.write('\x1b[1A')
# delete last line
sys.stdout.write('\x1b[2K')
print('imported data!')


#!! location of outputs
path = '../plots/'

pathexists = os.path.exists(path)
if pathexists:
    print(path + ' exists.')
else:
    print(path + ' does not exist.')
    os.makedirs(path)
    print(path + ' is created.')
plot_dir = os.path.abspath(sim_path+'/plot')
print('Output directory: ', plot_dir)


#!! comment out to plot them.
# # creating an StateVariables object
# img_StateVariables = statevariables.StateVariables(path, sim, ts, xyaver, yaver, ini=False, dyn=True)
# # plotting butterfly diagram
# img_StateVariables.plot(8,5)


# # creating an Rms object
# img_Rms = rms.Rms(path, sim, ts, xyaver, yaver, ini=False, dyn=True)
# # plotting butterfly diagram
# img_Rms.plot(8,5)


# creating an Butterfly object
img_Bfly = butterfly.Butterfly(path, sim, ts, xyaver, yaver, ini=False, dyn=True)
# plotting butterfly diagram
img_Bfly.plot()


# # creating an KOm object
# img_KOm = kom.KOm(path, 150, 670, 230, sim, ts, xyaver, yaver, dyn=True) #ini=False,
# # plotting k-omega diagram
# img_KOm.plot(detailed=True, cmap='hot_r')

#!! comment out the below portion to plot all the diagrams together.
# try:
#     try:
#         # creating an StateVariables object
#         img_StateVariables = statevariables.StateVariables(path, sim, ts, xyaver, yaver, ini=False, dyn=True)
#         # plotting butterfly diagram
#         img_StateVariables.plot(8,5)
#         print('saved statevariables')
#     except Exception as error:
#         # print("An exception occurred:", type(error))
#         print('\nAn exception of ', type(error), ' occured when plotting state variables!')
    
#     try:
#         # # creating an Rms object
#         img_Rms = rms.Rms(path, sim, ts, xyaver, yaver, ini=False, dyn=True)
#         # # plotting butterfly diagram
#         img_Rms.plot(8,5)
#         print('saved rms')
#     except Exception as error:
#         # print("An exception occurred:", type(error))
#         print('\nAn exception of ', type(error), ' occured when plotting rms values!')
    
#     try:
#         # # creating an Butterfly object
#         img_Bfly = butterfly.Butterfly(path, sim, ts, xyaver, yaver, ini=False, dyn=True)
#         # # plotting butterfly diagram
#         # img_Bfly.plot(ticks=[-.3, -.2, -.1, 0, .1, .2, .3])
#         img_Bfly.plot(ticks=[-.3, -.2, -.1, 0, .1, .2, .3])
#         print('saved butterfly')
#     except Exception as error:
#         # print("An exception occurred:", type(error))
#         print('\nAn exception of ', type(error), ' occured when plotting butterfly diagram!')
    
#     try:
#         # # creating an KOm object
#         img_KOm = kom.KOm(path, 150, 670, 230, sim, ts, xyaver, yaver, ini=False, dyn=True)
#         # # plotting k-omega diagram
#         img_KOm.plot()
#         print('saved kom')
#     except Exception as error:
#         # print("An exception occurred:", type(error))
#         print('\nAn exception of ', type(error), ' occured when plotting kom diagram!')
    

# finally:
#     print('\ncompleted.\n')
