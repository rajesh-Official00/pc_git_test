import phase
import modes
import numpy as np
import pencil as pc
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.convolution import convolve, Box1DKernel


path = './plots/'

sim = pc.get_sim(quiet=True)
ts = pc.read.ts()
xyaver = pc.read.aver(plane_list=['xy'])
yaver = pc.read.aver(plane_list='y')



kinematic = phase.Phase(path, 100, 2100, 225, 'ortho', sim, ts, xyaver, yaver, ini=False, dyn=True)
# saturated = phase.Phase(path, 4500, 7000, 225, 'ortho', sim, ts, xyaver, yaver, ini=False, dyn=True)


indx_kin = []
indx_sat = []


how_many = 3
hm = how_many
# for i in [2,3,4,5]:
for i in [1.5, 2, 2.5, 3]:
    indx_kin.append(kinematic.indx_k(i))
    # indx_sat.append(saturated.indx_k(i))

P_kin = kinematic.Power(indx_kin).T
# P_sat = saturated.Power(indx_sat).T

label_kin = []
P_kin_filt = np.empty(np.shape(P_kin))
# P_sat_filt = np.empty(np.shape(P_sat))
for i, indx in enumerate(indx_kin):
    P_kin_filt[i,:] =convolve(P_kin[i,:], Box1DKernel(5))
    label_kin.append(kinematic.k_til[:kinematic.upto_indx][indx])
    # P_sat_filt[i,:] =convolve(P_sat[i,:], Box1DKernel(5))


#! plot kinematic and saturated power vs kx
fig, axs = plt.subplots(2,1, sharex=True, figsize=(8,5))#sharex=True,

kinematic.plot(axs[0], P_kin[1], label='kinematic', c='k', alpha=0.8)
# kinematic.plot(axs[0], P_kin_filt[3], label='filtered')
# saturated.plot(axs[1], P_sat[1], label='saturated', c='k', alpha=0.8)
# saturated.plot(axs[1], P_sat_filt[3], label='filtered')

# axs[1].set_xlim(0,kinematic.om_til[kinematic.upto_indx-1])
axs[1].set_xlim(kinematic.xlim(kinematic.om_til))
# axs[1].set_ylim(0,1.75)
axs[0].grid()
axs[1].grid()
plt.savefig('test_plots/spectra.png')


f_om = kinematic.fmodes(kinematic.k_til[indx_kin], qq=True)

p0_om = kinematic.pmodes(kinematic.k_til[indx_kin], 0)
p1_om = kinematic.pmodes(kinematic.k_til[indx_kin], 1)



om_kin = kinematic.om_til[np.argmin(np.abs(kinematic.om_til-0)):]
# indx_f = np.argmin(np.abs(om_kin-f_om))
idl_f_kin = []
idu_f_kin = []
P_f_kin = []
om_f_kin = []

# d = [1.12, 1.5, 1.7, 2.0]
# u = [2.10, 2.2, 2.5, 2.6]

d = [0.8, 0.8, 1.0]
u = [1.6, 2.2, 2.6]

for i in range(hm):
    idl_f_kin.append(np.argmin(np.abs(om_kin - d[i])))
    idu_f_kin.append(np.argmin(np.abs(om_kin - u[i])))

    # P_f_kin.append(P_kin[i, idl_f_kin[i]:idu_f_kin[i]])
    P_f_kin.append(P_kin_filt[i, idl_f_kin[i]:idu_f_kin[i]])
    om_f_kin.append(om_kin[idl_f_kin[i]:idu_f_kin[i]])


fig, axs = plt.subplots(3,1, sharex=True, sharey=True, figsize=(8,6))#sharex=True,
# fig, axs = plt.subplots(4,1, sharey=True, figsize=(8,6))#sharex=True,

for i in range(3):
    kinematic.plot(axs[i], P_kin[i], c='k', alpha=0.4, label=fr'$\tilde{{k}}_x={round(label_kin[i],4)}$')
    axs[i].axvline(x=f_om[i], ls=':', c='b')
    axs[i].axvline(x=p0_om[i], ls=':', c='r')
    axs[i].axvline(x=p1_om[i], ls=':', c='r')

axs[1].set_xlim(kinematic.xlim(kinematic.om_til))
axs[0].set_ylim(0,)
axs[0].legend()

plt.xlabel(r"$\tilde{\omega}$")
fig.supylabel(r"$\tilde{P}(\tilde{\omega})$")
plt.suptitle(r"$kimnematic$")
#plt.yscale('log')
plt.tight_layout()
plt.savefig('test_plots/multi_k_kin.png')

fig, axs = plt.subplots(1,hm, sharey=True, figsize=(12,5))#sharex=True,

for i in range(hm):
    axs[i].plot(om_f_kin[i], P_f_kin[i], c='k', ls=':', label=fr'$\tilde{{k}}_x={round(label_kin[i],3)}$')
    axs[i].axvline(x=f_om[i], ls='-.', c='r', alpha=0.6)
    # axs[i].axvline(x=f_om[i], ls=':', c='b')
    # axs[i].set_xlim(om_f_kin[i][i], om_f_kin[i][-1])
    axs[i].legend()

axs[0].set_ylim(0,)

plt.xlabel(r"$\tilde{\omega}$")
fig.supylabel(r"$\tilde{P}(\tilde{\omega})$")
plt.suptitle(r"$kimnematic$")
#plt.yscale('log')
plt.tight_layout()
plt.savefig('test_plots/for_diff_kx.png')

def lorentzian(x, a, b, c, d, e):
    y = a/((x-b)**2+np.exp(c))+d+e*x
    return y

para_f_kin = np.zeros((hm,5))
para_f_sat = np.zeros((hm,5))

a_kin = [0.4, 0.4, 0.8]
# b_kin = [1.6, 1.7, 2.1, 2.3]
b_kin = [1.25, 1.55, 1.9]
d_kin = [-1.0, -1.0, -1.0]
e_kin = [0.1, 0.1, 0.1]

a_sat = [0.4, 0.4, 0.8]
# b_sat = [1.6, 1.7, 2.2, 2.3]
b_sat = [0.7, 0.9, 1.1]
d_sat = [-1.0, -1.0, -1.0]
e_sat = [0.1, 0.1, 0.1]


x_kin = np.zeros((hm,2000))
y_kin = np.zeros((hm,2000))
x_sat = np.zeros((hm,2000))
y_sat = np.zeros((hm,2000))

for i in range(hm):
    para_f_kin[i,:] = kinematic.mode_fit_para(lorentzian, om_f_kin[i], P_f_kin[i], p0 = np.array([a_kin[i],b_kin[i],-2,d_kin[i],e_kin[i]]))
    x_kin[i] = np.linspace(para_f_kin[i,1]-1.0, para_f_kin[i,1]+1.0, 2000)
    y_kin[i] = kinematic.mode_fit_extend(lorentzian, x_kin[i], para_f_kin[i,:], cont=True)

    # para_f_sat[i,:] = saturated.mode_fit_para(lorentzian, om_f_sat[i], P_f_sat[i], p0 = np.array([a_sat[i],b_sat[i],-2,d_sat[i],e_sat[i]]))
    # x_sat[i] = np.linspace(para_f_sat[i,1]-1.0, para_f_sat[i,1]+1.0, 2000)
    # y_sat[i] = saturated.mode_fit_extend(lorentzian, x_sat[i], para_f_sat[i,:])


fig, axes = plt.subplots(1,hm, sharey=True, figsize=(12,5))
for i in range(hm):
    axes[i].plot(om_f_kin[i], P_f_kin[i], color='k', ls=":", alpha=0.6)#, label=r'$kinematic$ $phase$')
    axes[i].plot(x_kin[i], y_kin[i], ls='--', color='k', label=fr'$\tilde{{k}}_x={round(label_kin[i],3)}$')#fitted_f_kin_lor-para_f_kin_lor[3],
    # axes[i].axvline(x=f_om[i], ls='-.', c='r', alpha=0.6)
    axes[i].legend(loc='lower left', bbox_to_anchor=(0.0, 1.0))
    axes[i].set_xlabel(r"$\tilde{\omega}$")

# axes[1].set_xlabel(r"$\tilde{\omega}$")
axes[0].set_ylabel(r"$\tilde{P}(\tilde{\omega})$")
plt.suptitle(r"$kinematic$")
plt.tight_layout()
plt.savefig('test_plots/for_diff_kx_fitted.png')


u_d_kin = 0.1794    #NOTE: the value is calculated by above method


mode_mass_kin = np.zeros(hm)
shift_kin = np.zeros(hm)
indx_fwhm_kin = np.zeros((hm,2))
ln_wd_kin = np.zeros(hm)

for i in range(hm):
    mode_mass_kin[i] = kinematic.mode_mass(x_kin[i], y_kin[i], u_d_kin)#; mode_mass_sat = mode_mass(f_sat[0], f_sat[1], u_d_sat)
    shift_kin[i] = kinematic.shift(para_f_kin[i,1], f_om[i])
    indx_fwhm_kin[i,:] = kinematic.indx_fwhm(y_kin[i], x_kin[i], para_f_kin[i,1])
    ln_wd_kin[i] = kinematic.line_width(x_kin[i,:], indx_fwhm_kin[i,:], f_om[i])

def st_line(x,a,b):
    return a*x+b

para_st_kin, _ = curve_fit(st_line, label_kin, shift_kin)

x_kin = np.linspace(label_kin[0], label_kin[-1], 100)
y_kin = st_line(x_kin, *para_st_kin)



fig, ax = plt.subplots(3,1, sharex=True, figsize=(6,5))
size = 10
ax[0].scatter(label_kin, mode_mass_kin, color='k', s=size)#fitted_f_kin_lor-para_f_kin_lor[3],
ax[1].scatter(label_kin, shift_kin, color='k', s=size)#fitted_f_kin_lor-para_f_kin_lor[3],
ax[2].scatter(label_kin, ln_wd_kin, color='k', s=size)#fitted_f_kin_lor-para_f_kin_lor[3],

ax[1].plot(x_kin, y_kin, c='r', ls=':')

ax[0].set_ylim(np.min(mode_mass_kin)-10, np.max(mode_mass_kin)+10)
ax[1].set_ylim(np.min(shift_kin)-0.02, np.max(shift_kin)+0.01)
ax[2].set_ylim(np.min(ln_wd_kin)-0.03, np.max(ln_wd_kin)+0.03)

ax[0].set_ylabel(r'$\mu_{f,kin}$')
ax[1].set_ylabel(r'$\delta\omega_{f,kin}^2/\omega_f^2$')
ax[2].set_ylabel(r'$\Gamma_{f,kin}$')
plt.xlabel(r'$\tilde{k}_x$')
plt.savefig('test_plots/modes.png', dpi=300)
# plt.show()