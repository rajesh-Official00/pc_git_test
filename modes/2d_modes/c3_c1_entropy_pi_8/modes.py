import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from peakutils import indexes

class modes():

    def __init__(self, sim, ts, xyaver, yaver, dyn=False, ini=False):
        self.sim = sim
        self.ts = ts
        self.xyaver = xyaver
        self.yaver = yaver

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
        #reading the latest thermodynamic variables
        self.rho = self.xyaver.xy.rhomz[-1,:]
        self.pre = self.xyaver.xy.ppmz[-1,:]
        self.tem = self.xyaver.xy.TTmz[-1,:]
        self.scale_height = self.cp*(1-1/self.gamma)*self.tem/self.gz
        if dyn==True:
            self.bxmz = self.xyaver.xy.bxmz[-1,:]
            self.bymz = self.xyaver.xy.bymz[-1,:]
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



    def index(x, y):
        return np.argmin(np.abs(x - y))

    def temperature(self, i, j):
        def const(t, a):
            return a
        z_ij = self.z[i:j]
        T_ij = self.tem[i:j]
    
        popt, _ = curve_fit(const, z_ij, T_ij)
        # z_fit = z[i:j]
        T_fit = const(z_ij, popt[0])
        return round(np.average(T_fit),3)
    
    def cs(self, gamma, R, T):
        return np.round(np.sqrt(gamma*R*T), 3)
    
    def urms(self, ti, tj):
        def g(t, a, b):
            return a*t+b
        popt, _ = curve_fit(g, self.ts.t[ti:tj], self.ts.urms[ti:tj])
        u_fit = g(self.ts.t[ti:tj], *popt)
        return round(np.average(u_fit),4)
    
    def D(self, L0, u_rms, cs_d):
        return np.round(L0*u_rms/cs_d, 3)
    
    def FT(self, uz_real, norm: str):
        return np.fft.fftn(uz_real, s=None, axes=(-2, -1), norm=norm)
    
    def logP(self, uz_fourier,D):
        return np.log(np.abs(uz_fourier/D**2))

    def omega_tilde(self, tt, t1, t2):
        "tt=yaver.t"
        t_gd = tt[t1:t2]  #time interval where urms has reached a steady state
        t_len = np.size(t_gd)
        dom = 2*np.pi/t_len #unit step alomg omega direction
        if t_len%2 == 0:
            fom = np.arange(0, t_len/2+1)
            rom = -np.flip(np.arange(1, t_len/2))
            om = np.concatenate((fom, rom))*dom
        else:
            fom = np.arange(0, t_len/2)
            rom = -np.flip(np.arange(1, t_len/2))
            om = np.concatenate((fom, rom))*dom
        return om/self.omega0
    
    def upto(self, tt, t1, t2):
        "tt=yaver.t"
        t_gd = tt[t1:t2]  #time interval where urms has reached a steady state
        t_len = np.size(t_gd)
        if t_len%2 == 0:
            fom = np.arange(0, t_len/2+1)
        else:
            fom = np.arange(0, t_len/2)
        return len(fom)

    def k_tilde(self):
        nx = len(self.grid.x)
        dkx = 2*np.pi/self.lx    #unit step along kx direction
        if nx%2 == 0:
            fnx = np.arange(0, nx/2+1)
            rnx = -np.flip(np.arange(1, nx/2))
            kx = np.concatenate((fnx, rnx))*dkx
        else:
            fnx = np.arange(0, nx/2)
            rnx = -np.flip(np.arange(1, nx/2))
            kx = np.concatenate((fnx, rnx))*dkx
        return kx*self.L0
    
    def fmodes(self,kx_tilde):
        # om_sq = self.gz*(kx_tilde/self.L0)*(1-self.q)/(1+self.q)
        om_sq = self.gz*(kx_tilde/self.L0)#*(1-self.q)/(1+self.q)
        return np.sqrt(om_sq/self.omega0**2)
    
    def pmodes(self,kx,n:int):
        om_sq = self.gz**2/(2*self.cs_d)**2+self.cs_d**2*((kx/self.L0)**2+((n+0.5)*np.pi/(9*self.lz/10))**2)
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
    
    def plot(x_data, y_data, **kwargs):
        plt.plot(x_data, y_data)
        plt.xlim(x_data[0], x_data[-1])

    pass