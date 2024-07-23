import matplotlib
import matplotlib.pyplot as plt
import modes
matplotlib.use('Agg')

# plt.rcParams.update({'font.size': 14})
# plt.rcParams['text.usetex'] = True

class StateVariables(modes.Modes):
    """
    Plot state variables as functions of z
    """

    def __init__(self, path, sim, ts, xyaver, yaver, ini=False, dyn=True):
        super().__init__(sim, ts, xyaver, yaver, ini=False, dyn=True)

        self.path = path

    def plot(self, width, height):
        """
        width: width of the figure
        height: height of the figure
        """
        
        fig, axs = plt.subplots(3,1, figsize=(width,height), sharex=True)
        axs[0].plot(self.z, self.rho, color='k')
        axs[1].plot(self.z, self.pre, color='k')
        axs[2].plot(self.z, self.tem, color='k')

        axs[0].set_yscale("log")
        axs[1].set_yscale("log")

        axs[2].set_xlim(self.z[0], self.z[-1])

        axs[0].set_ylabel(r"$\rho(z)$")
        axs[1].set_ylabel(r"$P(z)$")
        axs[2].set_ylabel(r"$T(z)$")

        # xticks = (np.pi/10)*np.array([-9, -7, -5, -3, -1, 0, 1])
        # axs[2].set_xticks(xticks)
        # axs[2].set_xticklabels([r'$-\frac{9\pi}{10}$', r'$-\frac{7\pi}{10}$', r'$-\frac{5\pi}{10}$', \
        #                      r'$-\frac{3\pi}{10}$', r'$-\frac{\pi}{10}$', r'$0$', r'$\frac{\pi}{10}$'])
        plt.xlabel(r"$z$")
        plt.tight_layout()
        plt.savefig(self.path+'statevariables.pdf',dpi=300)
    

    pass