import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import modes
matplotlib.use('Agg')

# plt.rcParams.update({'font.size': 14})
# plt.rcParams['text.usetex'] = True

class KOm(modes.Modes):
    """
    class to plot k-omega diagram.
    t1: initial time
    t2: final time
    z_ref: indx of z=0
    """
    def __init__(self, path, t1, t2, z_ref, sim, ts, xyaver, yaver, ini=False, dyn=True):
        super().__init__(sim, ts, xyaver, yaver, ini=False, dyn=True)

        self.path = path
        self.t1 = t1
        self.t2 = t2
        self.z_ref = z_ref

        self.indx_t1 = np.argmin(np.abs(self.yaver.t - self.t1))  #to avoid the transient effects below t=170
        self.indx_t2 = np.argmin(np.abs(self.yaver.t - self.t2))

        self.uz_real = self.yaver.y.uzmxz[self.indx_t1:self.indx_t2,:,self.z_ref]
        self.uz_fourier = super().FT(self.uz_real, 'ortho')
        self.log_P = super().logP(self.uz_fourier, self.d)
        self.om_til = super().omega_tilde(self.indx_t1, self.indx_t2)
        self.k_til = super().k_tilde()


    def plot(self, vmin=-3.5, cmap='afmhot_r', detailed=False):
        """
        function to plot k-omega diagram.
        """

        levels = np.linspace(-4, np.max(self.log_P), 1000)
        vmax = np.max(self.log_P)
    
        [X, Y] = np.meshgrid(self.k_til, self.om_til)

        fig, ax = plt.subplots(1, figsize=(7,5))
        # ax = fig.add_subplot(1,1,1)
        diag = ax.contourf(np.fft.fftshift(X), np.fft.fftshift(Y), np.fft.fftshift(self.log_P), levels=levels, cmap=cmap, vmin=vmin, vmax=vmax, extend='min')

        if detailed:
            plt.plot(self.xx, np.sqrt(self.gz*self.xx/(self.omega0**2*self.L0)), 'k')
            plt.plot(self.xx, np.sqrt(self.gz*self.xx/(self.omega0**2*self.L0)*(1-self.q)/(1+self.q)), ls='dotted', c='k')
            # plt.plot(self.k_til, np.sqrt(self.gz*self.k_til/(self.omega0**2*self.L0)), 'k')
            # plt.plot(self.k_til, np.sqrt(self.gz*self.k_til/(self.omega0**2*self.L0)*(1-self.q)/(1+self.q)), ls='dotted', c='k')
            plt.plot(self.k_til, self.cs_d*self.k_til/(self.omega0*self.L0), ls='--', c='k')
            plt.plot(self.k_til, self.cs_u*self.k_til/(self.omega0*self.L0), ls=':', c='k')

        plt.xlim(-2*np.pi,2*np.pi)
        cbar = fig.colorbar(diag, cax=None, ax=ax)
        ticks = np.arange(vmin, vmax)
        cbar.set_ticks(ticks)
        cbar.ax.set_xlabel(r'$lnP$', labelpad=20)
        plt.xlabel(r"$\tilde{k}_x$")
        plt.ylabel(r"$\tilde{\omega}$")
        plt.savefig(self.path+'k_om.png',dpi=300)
        # plt.show()
        return diag

    pass
