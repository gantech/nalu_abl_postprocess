import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import nalu_read_data as nrd
import abl_tools as at


class naluABLCase:
    Qo = 0.0 #Surface heating rate
    g = 9.81 #Acceleration due to gravity
    kappa = 0.41 #Kolmogorov constant
    To = 300.0 #Reference temperature
    nx = 126 #Number of points in the x - direction
    ny = 126 #Number of points in the y - direction
    zparts = ['zplane_0080.0', 'zplane_0320.0', 'zplane_0600.0'] #Parts corresponding
    nzSpectra = len(zparts) #Number of planes of data for spectra
    nzVar = 1 #Number of planes at which variances are computed. To be determined from reading mean/variance files
    Lx = 5000.0 #Length of domain in x-direction
    Ly = 5000.0 #Length of domain in y-direction
    dx = Lx/nx #Resolution in x-direction
    dy = Ly/ny #Resolution in y-direction
    nt = 7     #Number of time-steps for analysis
    nWavenumbers = np.int(nx/np.sqrt(2))
    
    
    zvar = np.empty(1)
    umean = np.empty((1,nt))
    vmean = np.empty((1,nt))
    wmean = np.empty((1,nt))
    thetamean = np.empty((1,nt))
    uvar_res = np.empty((1,nt))
    vvar_res = np.empty((1,nt))
    wvar_res = np.empty((1,nt))
    upvp_res = np.empty((1,nt))
    upwp_res = np.empty((1,nt))
    vpwp_res = np.empty((1,nt))
    uvar_sfs = np.empty((1,nt))
    vvar_sfs = np.empty((1,nt))
    wvar_sfs = np.empty((1,nt))
    upvp_sfs = np.empty((1,nt))
    upwp_sfs = np.empty((1,nt))
    vpwp_sfs = np.empty((1,nt))
    wpthetap_res = np.empty((1,nt))
    wpthetap_sfs = np.empty((1,nt))
    wpcube_res = np.empty((1,nt))
    
    uspectra = np.zeros((nWavenumbers, nt, nzSpectra))
    vspectra = np.zeros((nWavenumbers, nt, nzSpectra))
    wspectra = np.zeros((nWavenumbers, nt, nzSpectra))
    thetaspectra = np.zeros((nWavenumbers, nt, nzSpectra))
    
    utau = np.empty(nt)
    L = np.empty(nt)

    def __init__(self):
        pass
    
    def process_zvar(self):
        self.zvar = nrd.get_heights_var('abl_T_stats.dat')
        self.nzVar = np.size(self.zvar)
        
    def process_theta_mean(self):
        self.thetamean.resize((self.nzVar,self.nt))
        self.thetamean = nrd.read_theta_mean('abl_T_stats.dat', self.nzVar)
        
    def process_vel_mean(self):
        self.umean.resize((self.nzVar,self.nt))
        self.vmean.resize((self.nzVar,self.nt))
        self.wmean.resize((self.nzVar,self.nt))
        self.umean, self.vmean, self.wmean = nrd.read_vel_mean('abl_U_stats.dat', self.nzVar) 
        
    def process_variances(self):
        self.uvar_res.resize((self.nzVar,self.nt))
        self.vvar_res.resize((self.nzVar,self.nt))
        self.wvar_res.resize((self.nzVar,self.nt))
        self.upvp_res.resize((self.nzVar,self.nt))
        self.upwp_res.resize((self.nzVar,self.nt))
        self.vpwp_res.resize((self.nzVar,self.nt))
        self.uvar_sfs.resize((self.nzVar,self.nt))
        self.vvar_sfs.resize((self.nzVar,self.nt))
        self.wvar_sfs.resize((self.nzVar,self.nt))
        self.upvp_sfs.resize((self.nzVar,self.nt))
        self.upwp_sfs.resize((self.nzVar,self.nt))
        self.vpwp_sfs.resize((self.nzVar,self.nt))
        self.wpthetap_res.resize((self.nzVar,self.nt))
        self.wpthetap_sfs.resize((self.nzVar,self.nt))
        self.wpcube_res.resize((self.nzVar,self.nt))

        #    umean, vmean, wmean, thetamean, uvar_res, vvar_res, wvar_res, upvp_res, upwp_res, vpwp_res, uvar_sfs, vvar_sfs, wvar_sfs, upvp_sfs, upwp_sfs, vpwp_sfs, wpthetap_res, wpthetap_sfs, wpcube_res = nrd.read_mean_var()
        self.uvar_res, self.vvar_res, self.wvar_res, self.upvp_res, self.upwp_res, self.vpwp_res, self.thetavar_res, self.wpthetap_res, self.wpcube_res = nrd.read_var('abl_var_stats.dat', self.nzVar)
    
    def process_spectra(self):
        ds = nc.Dataset('zplaneData.e')
        t = np.zeros(self.nt)
        for iz,zplane in enumerate(self.zparts):
            for tidx in range(self.nt):
                ux, uy, uz, tN = nrd.read_nalu_plane_vel(self.nx,self.ny,tidx,ds,zplane)
                self.uspectra[:,tidx,iz] = at.calc_twod_spectra(ux)
                self.vspectra[:,tidx,iz] = at.calc_twod_spectra(uy)
                self.wspectra[:,tidx,iz] = at.calc_twod_spectra(uz)
#                theta, tN = nrd.read_nalu_plane_theta(self.nx,self.ny,tidx,ds,zplane) 
#                self.thetaspectra[:,tidx,iz] = at.calc_twod_spectra(theta)
                
    def process_utau(self):
        utau = nrd.read_utau('abl_uTau_stats.dat')
    
    def process_zi(self):
        self.zi = at.calc_zi(self.wpthetap_res, self.zvar)
    
    def process_L(self):
        self.L.resize(self.nt)
        if (self.Qo > 0):
            self.L = at.calc_L(self.utau, self.g, self.To, self.Qo, self.kappa)
        else:
            self.L[:] = np.inf
    

    def plot_mean_profile(self, phimean, z, axisLabel):
        fig = plt.figure()
        plt.plot(phimean, z)
        plt.xlabel(axisLabel)
        plt.ylabel('z')
        plt.savefig('{}_profile.png'.format(axisLabel))
            
    def plot_latest_profiles(self):
        self.plot_mean_profile(self.umean[-1,:], self.zvar, 'u')
        self.plot_mean_profile(self.vmean[-1,:], self.zvar, 'v')    
        self.plot_mean_profile(self.vmean[-1,:], self.zvar, 'w')
        self.plot_mean_profile(self.thetamean[-1,:], self.zvar, 'theta')
        self.plot_mean_profile(self.uvar_res[-1,:], self.zvar, 'uvar_res')
        self.plot_mean_profile(self.vvar_res[-1,:], self.zvar, 'vvar_res')
        self.plot_mean_profile(self.wvar_res[-1,:], self.zvar, 'wvar_res')
        self.plot_mean_profile(self.upvp_res[-1,:], self.zvar, 'upvp_res')
        self.plot_mean_profile(self.upwp_res[-1,:], self.zvar, 'upwp_res')
        self.plot_mean_profile(self.vpwp_res[-1,:], self.zvar, 'vpwp_res')
        self.plot_mean_profile(self.wpthetap_res[-1,:], self.zvar, 'wpthetap_res')
        self.plot_mean_profile(self.wpcube_res[-1,:], self.zvar, 'wpcube_res')



    def plot_single_field_spectra(self, kr, phiSpectra, axisLabel, zplane):
        fig = plt.figure()
        plt.plot(kr, phiSpectra)
        plt.xlabel('Wavenumber k')
        plt.ylabel('E - {}'.format(axisLabel))
        plt.savefig('{}_{}_spectra.png'.format(axisLabel,zplane))
        
        
    def plot_latest_spectra(self):
        kr = np.arange(self.nWavenumbers)*2*np.pi/self.Lx #Will not work when Lx != Ly
        for iz, zplane in enumerate(self.zparts):
            self.plot_single_field_spectra(kr, self.uspectra[:,-1,iz], 'Ux', zplane)
            self.plot_single_field_spectra(kr, self.vspectra[:,-1,iz], 'Uy', zplane)
            self.plot_single_field_spectra(kr, self.wspectra[:,-1,iz], 'Uz', zplane)
        
        
if __name__=="__main__":

    testCase = naluABLCase()
    testCase.process_zvar()
    testCase.process_theta_mean()
    testCase.process_vel_mean()
    testCase.process_variances()
    testCase.process_spectra()
    testCase.process_utau()
    testCase.process_zi()
    testCase.process_L()
    testCase.plot_latest_profiles()
    testCase.plot_latest_spectra()

