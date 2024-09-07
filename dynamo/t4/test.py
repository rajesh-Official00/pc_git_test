# %%
import pencil as pc
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# %%
ts=pc.read.ts()

# %%
#for exp fitting
def func(t, a, b):
    return a*np.exp(t*b)

i = np.argmin(abs(ts.t-200))
j = np.argmin(abs(ts.t-900))

t_exp = ts.t[i:j]
b_exp = ts.brms[i:j]

popt, pcov = curve_fit(func, t_exp, b_exp, [0,0.0016])

b_fit = func(ts.t[i:j], popt[0], popt[1])

# %%
# popt

# %%
plt.figure(figsize=(10,5))
plt.plot(ts.t, ts.brms)
plt.plot(t_exp, b_fit, linewidth=2.5, linestyle = '--')
plt.yscale("log")
# plt.savefig('output.jpg')
plt.grid()
plt.show()

# # %%
# (10**(-2)-10**(-4))/(375)


# # %%
# x = np.linspace(0,10)
# plt.plot(x, np.exp(x))
# plt.yscale("log")


# # %%



