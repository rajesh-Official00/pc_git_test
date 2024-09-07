import matplotlib.pyplot as plt
import numpy as np
import pencil as pc
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit
from scipy.integrate import trapezoid as trap
from astropy.convolution import convolve, Box1DKernel
import modes
# from IPython.display import display, Latex

# plt.rcParams.update({'font.size': 11})
# plt.rcParams['text.usetex'] = True

# sim = pc.get_sim(quiet=True)
# ts = pc.read.ts()

# xyaver = pc.read.aver(plane_list=['xy'])

# yaver = pc.read.aver(plane_list='y')


class Phase(modes.modes):
    """
    To analyze different phases
    """

    def __init__(self, path, t1, t2, z_ref, norm, sim, ts, xyaver, yaver, ini=False, dyn=True):
        super().__init__(sim, ts, xyaver, yaver, ini=False, dyn=True)

        self.path = path
        self.t = self.yaver.t
        self.z_ref = z_ref
        self.norm = norm

        self.kf = super().get_kf()
        self.indx_t1 = np.argmin(np.abs(self.yaver.t-t1))
        self.indx_t2 = np.argmin(np.abs(self.yaver.t-t2))

        self.uz_real = self.yaver.y.uzmxz[self.indx_t1:self.indx_t2,:,self.z_ref]
        self.uz_fourier = super().FT(self.uz_real, self.norm)
        # self.uz_fourier = super().FT(self.uz_real, 'ortho')
        self.log_P = super().logP(self.uz_fourier, self.d)
        self.om_til = super().omega_tilde(self.indx_t1, self.indx_t2)
        self.k_til = super().k_tilde()
        self.upto_indx = super().upto(self.indx_t1, self.indx_t2)

    def indx_k(self, k_xtil:np.ndarray):
        indx = np.argmin(np.abs(self.k_til-k_xtil))
        return indx
    
    def Power(self, indx:int):
        # P = np.exp(self.log_P[:self.upto_indx,indx])
        P = np.exp(self.log_P[np.argmin(np.abs(self.om_til-0)):,indx])
        return P
    
    def Power_filtered(self, power:np.ndarray, sigma):
        # P_filt = gaussian_filter(power[:self.upto_indx],sigma)
        P_filt = gaussian_filter(power[np.argmin(np.abs(self.om_til-0)):],sigma)
        return P_filt
    
    def f_freq(self, k_tilx, qq=False):
        # freq = round(self.fmodes(k_tilx),3)
        freq = self.fmodes(k_tilx, qq)
        return freq
    
    def p_freq(self, k_tilx, num):
        if isinstance(num, int):
            freq = self.pmodes(k_tilx,num)
        if isinstance(num, list):
            freq=[]
            for i in range(len(num)):
                freq.append(self.pmodes(k_tilx,i))
        return freq
    
    def plot(self, ax, P:np.ndarray, **kwargs):
        ax.plot(self.om_til[np.argmin(np.abs(self.om_til-0)):], P, **kwargs)
        # ax.grid()
        # ax.set_ylim(0,)
        ax.legend()

    def xlim(self, x_data:np.ndarray):
        x_data = x_data[np.argmin(np.abs(self.om_til-0)):]
        min = x_data[0]
        max = x_data[-1]
        return (min,max)
    
    def cs_du(self, k_tilx):
        csd = self.cs_d*k_tilx/(self.L0*self.omega0)
        csu = self.cs_u*k_tilx/(self.L0*self.omega0)
        return (csd, csu)
    
    def mode_mass(self, x_data, y_data, u_d, **kwargs):
        """
        Calculate mode mass by integrating
        y_data in the interval x_data
        x_data: \omega_tilde range for a mode
        y_data: corresponding P(\omega_tilde)
        u_d: rms velocity in lower layer
        """
        # if kinematic == True:
        #     u_d = 0.1244
        # elif saturated == True:
        #     u_d = 0.0988
        # else:
        #     u_d = 1     #NOTE: not normalized
        norm = u_d/(3*self.kf)
        mass = trap(y_data, x=x_data)
        norm_mass = mass/norm
        return np.round(norm_mass, 4)

    def shift(self, f_om_calc, f_om):
        """
        Calculate relative freq shift from
        calculated \omega_f and theoretical
        \omega_f
        """
        return np.round((f_om_calc**2-f_om**2)/f_om**2, 4)

    def indx_fwhm(self, P_om, om_f, peak_f):
        """
        P_om: fitted spectrum (i.e., P(\omega)
        of P(\omega) vs \omega )
        om_f: freq range around \omega_f
        peak_f: freq of f mode calculated
        from fitting
        """
        indx_peak = np.argmin(np.abs(om_f-peak_f))
        peak = P_om[indx_peak]
        lh = np.argmin(np.abs(P_om[:indx_peak]-peak/2))
        uh = indx_peak + np.argmin(np.abs(P_om[indx_peak+1:]-peak/2))
        return [lh, uh]
    
    def line_width(self, om_f, indx_fwhm, f_om):
        return np.round((om_f[int(indx_fwhm[1])] - om_f[int(indx_fwhm[0])])/f_om, 4)

    def mode_fit_para(self, func, x_data, y_data, **kwargs):
        para, _ = curve_fit(func, x_data, y_data, **kwargs)
        return para
    
    def mode_fit_extend(self, func, x, para, cont=True):
        if cont:
            y = func(x, *para)
        else:
            y = func(x, *para)-x*para[4]-para[3]
        return y

    pass