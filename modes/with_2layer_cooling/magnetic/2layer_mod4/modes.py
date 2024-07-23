import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from peakutils import indexes

class modes():

    def __init__(self, sim, ts, xyaver, yaver, dyn=False, ini=True):
        self.sim = sim
        self.ts = ts
        self.xyaver = xyaver
        self.yaver = yaver

        self.nxgrid = sim.get_value('nxgrid')
        self.nygrid = sim.get_value('nygrid')
        self.nzgrid = sim.get_value('nzgrid')

        self.grid = sim.grid
        self.param = sim.param
        self.z = self.grid.z
        # index of z=0
        self.indx_zref = np.argmin(abs(self.z-0))
        # box geometry
        self.lx, self.ly, self.lz = self.param.get('lxyz')
        # other parameters
        self.gz = -self.param.get('gravz')
        self.cp = self.param.get('cp')
        self.gamma = np.round(self.param.get('gamma'), 3)
        self.R = np.round(self.cp*(1-(1/self.gamma)), 4)
        self.nu = self.param.get('nu')
        #reading the latest thermodynamic variables
        self.rho = self.xyaver.xy.rhomz[-1,:]
        self.pre = self.xyaver.xy.ppmz[-1,:]
        self.tem = self.xyaver.xy.TTmz[-1,:]
        self.scale_height = self.cp*(1-1/self.gamma)*self.tem/self.gz
        if dyn==True:
            self.bxmz = self.xyaver.xy.bxmz[-1,:]
            self.bymz = self.xyaver.xy.bymz[-1,:]
            try:
                self.bx2mz = self.xyaver.xy.bx2mz[-1,:]
                self.by2mz = self.xyaver.xy.by2mz[-1,:]
            except AttributeError:
                print("Have you included bx2mz and by2mz in xyaver.in?")

        if ini==False:
            # finding T_d
            self.T_d = self.temperature(10,self.indx_zref-10)
            # finding T_u
            self.T_u = self.temperature(self.indx_zref+10,-4)
            self.q = self.T_d/self.T_u
            # finding u_rms
            self.u_rms = self.urms(200, 800)
            #sound speed in upper and lower layer
            self.cs_d = self.cs(self.gamma, self.R, self.T_d)
            self.cs_u = self.cs(self.gamma, self.R, self.T_u)
            
            self.L0 = np.round(self.cs_d**2/self.gz, 3)
            self.omega0 = np.round(self.gz/self.cs_d, 3)
            self.d = self.D(self.L0,self.u_rms,self.cs_d)



    def get_kf(self):
        """
        This function will work for k.dat files
        created by kdat.py
        """
        f = open('./k.dat', 'r')
        first_line = f.readline()
        # print(first_line)
        kf = first_line.split(' ')[-1].strip()
        # print(type(kf))
        return float(kf)

    def u_rms_d(self, var):#, prnt=False):
        uu = var.uu
        if uu[1].shape == (self.nzgrid,self.nygrid,self.nxgrid):
            uu_d = var.uu[:,:self.indx_zref,:,:]
        else:
            uu_d = var.uu[:,3:self.indx_zref+3,3:-3,3:-3]
        # if prnt: print(uu_d.shape)
        uu_d_sq = uu_d**2
        uu_d_sq_sum = np.sum(uu_d_sq, axis=0)
        u_rms_d = np.sqrt(np.average(uu_d_sq_sum))
        u_rms_d = round(u_rms_d, 4)
        return u_rms_d

    def reynolds(self, var):
        """
        Calculate the reynoylds number of lower
        layer from var.dat file
        """
        Re = round(self.u_rms_d(var)/(self.nu*self.get_kf()), 4)
        return Re
    
    def ave_reynolds(self, *vars):
        """
        Calculate the average Reynoylds number of lower
        layer from the given multiple or single VAR files
        """
        rey = np.empty(len(vars))
        for i in range(len(vars)):
            rey[i] = self.reynolds(vars[i])
        rey = np.average(rey)
        rey = round(rey,4)
        return rey

    def mach(self, var):
        """
        Calculate the Mach number of lower
        layer from var.dat file
        """
        # print(self.u_rms_d(var))
        M = self.u_rms_d(var)/self.cs_d
        M = round(M, 4)
        return M
    
    def index(x, y):
        return np.argmin(np.abs(x - y))

    def temperature(self, i: int, j: int):
        """
        Calculate the average temperature of a layer
        from xyaverage.dat
        """
        def const(t, a):
            return a
        z_ij = self.z[i:j]
        T_ij = self.tem[i:j]
    
        popt, _ = curve_fit(const, z_ij, T_ij)
        # z_fit = z[i:j]
        T_fit = const(z_ij, popt[0])
        return round(np.average(T_fit),3)
    
    def cs(self, gamma, R, T):
        """
        Calculate the sound speed of a layer
        """
        return np.round(np.sqrt(gamma*R*T), 3)
    
    def urms(self, ti: int, tj: int):
        """
        Calculate the average temperature of a layer
        from time_series.dat
        """
        def g(t, a, b):
            return a*t+b
        popt, _ = curve_fit(g, self.ts.t[ti:tj], self.ts.urms[ti:tj])
        u_fit = g(self.ts.t[ti:tj], *popt)
        return round(np.average(u_fit),4)
    
    def D(self, L0, u_rms, cs_d):
        # return np.round(L0*u_rms/cs_d, 3)
        return np.round(L0*0.1245/cs_d, 3) # FIXME: u_rms =0.1205
    
    def FT(self, uz_real, norm: str):
        """
        Calculate the fourier transform of u_z along the last
        two axis i.e., axes=(-2,-1)
        """
        return np.fft.fftn(uz_real, s=None, axes=(-2, -1), norm=norm)
    
    def logP(self, uz_fourier,D):
        """
        Calculate the spectral power from FT(u_z)
        """
        return np.log(np.abs(uz_fourier/D**2))

    def omega_tilde(self, t1, t2):
        """
        Calculate the omega tilde from yaver.t
        """
        t_gd = self.yaver.t[t1:t2]  #time interval where urms has reached a steady state
        t_len = np.size(t_gd)
        dom = 2*np.pi/t_len #unit step alomg omega direction
        if t_len%2 == 0:
            fom = np.arange(0, t_len/2)
            rom = -np.flip(np.arange(1, t_len/2+1))
            # om = np.concatenate((fom, rom))*dom
            om = np.concatenate((rom, fom))*dom
        else:
            fom = np.arange(0, t_len/2)
            rom = -np.flip(np.arange(1, t_len/2))
            # om = np.concatenate((fom, rom))*dom
            om = np.concatenate((rom, fom))*dom
        return om/self.omega0
        # return om
    
    def upto(self, t1, t2):
        """
        Calculate the index of yaver.t for discarding half of the data.
        ref: Nyquist's theorem.
        """
        "tt=yaver.t"
        t_gd = self.yaver.t[t1:t2]  #time interval where urms has reached a steady state
        t_len = np.size(t_gd)
        if t_len%2 == 0:
            # fom = np.arange(0, t_len/2+1)
            fom = np.arange(0, t_len/2)
        else:
            fom = np.arange(0, t_len/2)
        return len(fom)

    def k_tilde(self):
        """
        Calculate the k_x tilde from yaver.t
        """
        nx = len(self.grid.x)
        dkx = 2*np.pi/self.lx    #unit step along kx direction
        if nx%2 == 0:
            fnx = np.arange(0, nx/2)
            rnx = -np.flip(np.arange(1, nx/2+1))
            # kx = np.concatenate((fnx, rnx))*dkx
            kx = np.concatenate((rnx, fnx))*dkx
        else:
            fnx = np.arange(0, nx/2)
            rnx = -np.flip(np.arange(1, nx/2))
            # kx = np.concatenate((fnx, rnx))*dkx
            kx = np.concatenate((rnx, fnx))*dkx
        # return kx*self.L0
        return kx
    
    def fmodes(self, kx_tilde):
        """
        Calculate the freq of f-mode.
        """
        # om_sq = self.gz*(kx_tilde/self.L0)*(1-self.q)/(1+self.q)
        om_sq = self.gz*(kx_tilde/self.L0)#*(1-self.q)/(1+self.q)
        return np.sqrt(om_sq/self.omega0**2)
    
    def pmodes(self, kx_tilde, n: int):
        """
        Calculate the freq of p-mode.
        n: order of p-modes, 0, 1, 2, 3, ...
        """
        om_sq = self.gz**2/(2*self.cs_d)**2+self.cs_d**2*((kx_tilde/self.L0)**2+((n+0.5)*np.pi/(9*self.lz/10))**2)
        return np.sqrt(om_sq/self.omega0**2)
    
    def pmodes_test(self,kx,n:int):
        om_sq = self.gz**2/(2*self.cs_d)**2+self.cs_d**2*((kx/self.L0)**2+((n+0.5)*np.pi/self.lz)**2)
        return np.sqrt(om_sq/self.omega0**2)
    
    def mode_finder(self, x_data, y_data, plot_it=False, thres=0.3, min_dist=100):
        indx = indexes(y_data, thres=thres, min_dist=min_dist)
        if plot_it:
            plt.figure(figsize=(14,5))
            plt.plot(x_data, y_data)
            plt.scatter(x_data[indx], y_data[indx], c='r')
            plt.xlim(x_data[0], x_data[-1])
            plt.ylim(0,)
            plt.show()
        else:
            return indx
        
    def mode_indx(self, index: int, lower_bound: int, upper_bound: int):
        lb = index-lower_bound
        ub = index+upper_bound
        indx_dict = {'lb': lb, 'ub': ub}
        return indx_dict
    
    def mode_data(self, y_data, x_data, mode_indx: dict):
        lb = mode_indx['lb']
        ub = mode_indx['ub']
        data_dict = {'x': x_data[lb:ub], 'y': y_data[lb:ub]}
        return data_dict

    def mode_data1(self, y_data, x_data, index: int, lower_bound: int, upper_bound: int):
        lb = index-lower_bound
        ub = index+upper_bound
        data_dict = {'x': x_data[lb:ub], 'y': y_data[lb:ub]}
        return data_dict
    
    def sigma(self, power, filtered_power):
        sig_tot = (power-filtered_power)**2
        sig_ave = np.sqrt(sum(sig_tot)/len(sig_tot))
        # sig = sig_ave*np.ones(len(sig_tot_kin[idl_f_kin:idu_f_kin]))
        return sig_ave

    def mode_fit(self, func, x_data, y_data, base=True, **kwargs):
        para, _ = curve_fit(func, x_data, y_data, **kwargs)
        y_fit = func(x_data, *para)
        if base:
            y = y_fit
        else:
            y = y_fit-x_data*para[4]-para[3]
        return y
    
    def plotter(ax, x_data, y_data, **kwargs):
        """
        Plot x_data vs y_data on an axis
        """
        a = ax.plot(x_data, y_data, **kwargs)
        ax.set_xlim(x_data[0], x_data[-1])
        return a
    
    def butterfly(self, fig, axs, save=False):
        self.bxmz = self.xyaver.xy.bxmz
        self.bymz = self.xyaver.xy.bymz
        t = self.xyaver.xy.t
        by = [self.bxmz, self.bymz]
        try:
            self.bx2mz = self.xyaver.xy.bx2mz
            self.by2mz = self.xyaver.xy.by2mz
            by2 = [self.bx2mz, self.by2mz]
            by2sqr = np.sqrt(by2)

            [X, Y] = np.meshgrid(t, self.grid.z)
            for i in range(len(axs.flat)):
                ims = axs[i].contourf(X, Y, np.transpose(by[i]/by2sqr[i]), 100, cmap='plasma', vmin=np.min(by/by2sqr), vmax=np.max(by/by2sqr))

            axs[0].set_xlabel(r'$t$')
            axs[0].set_ylabel(r'$z$')
            axs[0].set_title(r'$\langle B_x \rangle _{xy}/\langle B_y^2 \rangle _{xy}^\frac{1}{2}$')
            axs[0].set_xlim(0,)

            axs[1].set_xlabel(r'$t$')
            axs[1].set_ylabel(r'$z$')
            axs[1].set_title(r'$\langle B_y \rangle _{xy}/\langle B_y^2 \rangle _{xy}^\frac{1}{2}$')
            axs[1].set_xlim(0,)

            #cbar_ax = fig.add_axes([1.05, 0.15, 0.05, 0.7])
            fig.colorbar(ims, ax=axs.ravel().tolist())
            #     # plt.tight_layout()
            if save==True:
                plt.savefig('ave.jpg')
        except AttributeError:
            print("Not normalized!")

            [X, Y] = np.meshgrid(t, self.grid.z)
            for i in range(len(axs.flat)):
                ims = axs[i].contourf(X, Y, np.transpose(by[i]), 100, cmap='plasma', vmin=np.min(by/by2sqr), vmax=np.max(by/by2sqr))

            axs[0].set_xlabel(r'$t$')
            axs[0].set_ylabel(r'$z$')
            axs[0].set_title(r'$\langle B_x \rangle _{xy}$')
            axs[0].set_xlim(0,)

            axs[1].set_xlabel(r'$t$')
            axs[1].set_ylabel(r'$z$')
            axs[1].set_title(r'$\langle B_y \rangle _{xy}$')
            axs[1].set_xlim(0,)

            #cbar_ax = fig.add_axes([1.05, 0.15, 0.05, 0.7])
            fig.colorbar(ims, ax=axs.ravel().tolist())
            #     # plt.tight_layout()
            if save==True:
                plt.savefig('ave.jpg')
        
        

    pass