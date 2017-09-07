import netCDF4 as nc
import numpy as np

def get_velocity(ds, vNames, nsNames, ns, it):
    """ Get the velocity from dataset 'ds' on the nodeset 'ns' with the variable names 'vNames' and nodeset names 'nsNames' at time index 'it' """

    nsIn = nsNames.index("{}".format(ns))

    uxIn = vNames.index("velocity_x")
    uyIn = vNames.index("velocity_y")
    uzIn = vNames.index("velocity_z")

    ux = ds.variables['vals_nset_var{0:d}ns{1:d}'.format(uxIn+1, nsIn+1)][it]
    uy = ds.variables['vals_nset_var{0:d}ns{1:d}'.format(uyIn+1, nsIn+1)][it]
    uz = ds.variables['vals_nset_var{0:d}ns{1:d}'.format(uzIn+1, nsIn+1)][it]
    return [ux, uy, uz]

def get_temperature(ds, vNames, nsNames, ns, it):
    """ Get the temperature from dataset 'ds' on the nodeset 'ns' with the variable names 'vNames' and nodeset names 'nsNames' at time index 'it' """

    nsIn = nsNames.index("{}".format(ns))

    thetaIn = vNames.index("temperature")

    theta = ds.variables['vals_nset_var{0:d}ns{1:d}'.format(thetaIn+1, nsIn+1)][it]
    return  theta

def get_coordinates(ds, nsNames, ns):
    """ Get the coordinates from dataset 'ds' for the nodeset 'ns' """
    nsIn = nsNames.index(ns)
    nsNodes = ds.variables['node_ns{0:d}'.format(nsIn+1)][:] - 1
    x = ds.variables['coordx'][nsNodes]
    y = ds.variables['coordy'][nsNodes]
    z = ds.variables['coordz'][nsNodes]
    return [x, y, z]

def read_nalu_plane_vel(nx,ny,tidx,ds,ns):
    nsNames = ["%s" % nc.chartostring(nn)
          for nn in ds.variables['ns_names'][:]]
    vn = ["%s" % nc.chartostring(nn)
          for nn in ds.variables['name_nset_var'][:]]
    timeSteps = ds.variables['time_whole'][:]
    uxC, uyC, uzC = get_velocity(ds, vn, nsNames, ns, tidx)
    uxC = uxC.reshape(ny,nx)[:,:]
    uyC = uyC.reshape(ny,nx)[:,:]
    uzC = uzC.reshape(ny,nx)[:,:]
    
    time = timeSteps[tidx]
    
    return [uxC, uyC, uzC, time]

def read_nalu_plane_theta(nx,ny,tidx,ds,ns):
    nsNames = ["%s" % nc.chartostring(nn)
          for nn in ds.variables['ns_names'][:]]
    vn = ["%s" % nc.chartostring(nn)
          for nn in ds.variables['name_nset_var'][:]]
    timeSteps = ds.variables['time_whole'][:]
    theta = get_temperature(ds, vn, nsNames, ns, tidx)
    theta = thetaC.reshape(ny,nx)[:,:]
    
    time = timeSteps[tidx]
    
    return [theta, time]

def get_heights_var(filename):
    """Read the file containing mean profile of only quantity and determine the heights at which they are computed"""
    heights = np.genfromtxt(filename, comments='*', delimiter=',', max_rows=1)[1:-1]
    return heights

def read_theta_mean(filename, nheights):
    """Read the files containing mean temperature profile """
    tstats = np.loadtxt(filename, delimiter=',', usecols=range(nheights+1), comments='#')
    tmean = tstats[:,1:]
    return tmean

def read_vel_mean(filename, nheights):
    """Read the files containing mean velocity profiles """
    ustats = np.loadtxt(filename, delimiter=',', usecols=range(nheights*3+1), comments='#')
    umean = ustats[:,1::3]
    vmean = ustats[:,2::3]
    wmean = ustats[:,3::3]
    return [umean, vmean, wmean]

def read_var(filename, nheights):
    """Read the files containing variance profiles """
    varstats = np.loadtxt(filename, delimiter=',', usecols=range(nheights*9+1), comments='#')
    uvar_res = varstats[:,1::9]
    vvar_res = varstats[:,2::9]
    wvar_res = varstats[:,3::9]
    upvp_res = varstats[:,4::9]    
    upwp_res = varstats[:,5::9]
    vpwp_res = varstats[:,6::9]
    thetavar_res = varstats[:,7::9]
    wpthetap_res = varstats[:,8::9]
    wpcube_res = varstats[:,9::9]
    return [uvar_res, vvar_res, wvar_res, upvp_res, upwp_res, vpwp_res, thetavar_res, wpthetap_res, wpcube_res]
#    return [uvar_res, vvar_res, wvar_res, upvp_res, upwp_res, vpwp_res, uvar_sfs, vvar_sfs, wvar_sfs, upvp_sfs, upwp_sfs, vpwp_sfs, wpthetap_res, wpthetap_sfs, wpcube_res]
    
def read_utau(filename):
    """Read the file containing the history of utau """
    utau = np.loadtxt(filename, delimiter=',', usecols=range(2), comments='#')
    return utau
