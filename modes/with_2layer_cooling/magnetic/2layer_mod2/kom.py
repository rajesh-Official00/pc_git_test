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


    def plot(self):
        """
        function to plot k-omega diagram.
        """

        [X, Y] = np.meshgrid(self.k_til, self.om_til)
        diag = plt.contourf(np.fft.fftshift(X), np.fft.fftshift(Y), np.fft.fftshift(self.log_P), 4000, cmap='afmhot_r', vmin=-4, vmax=np.max(self.log_P))
        plt.xlim(-2*np.pi,2*np.pi)
        plt.colorbar()
        plt.xlabel(r"$\tilde{k}_x$")
        plt.ylabel(r"$\tilde{\omega}$")
        plt.savefig(self.path+'k_om.png',dpi=300)
        # plt.show()
        return diag

    pass