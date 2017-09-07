import numpy as np

def get_mean_fluct(u):
    ny,nx = np.shape(u)
    umean = np.zeros((ny,nx))
    umean[:,:] = np.mean(u[:,:])
    ufluct = u - umean
    return [umean, ufluct]

def get_twod_sp_count(nx, ny, nWavenumbers):
    spCount = np.empty(nWavenumbers)
    k1 = np.fft.fftfreq(nx)*nx
    k2 = np.fft.fftfreq(ny)*ny
    for i in range(nx):
        for j in range(ny):
            r = np.rint(np.sqrt(k1[i]*k1[i] + k2[j]*k2[j]))
            if(r < nWavenumbers):
                spCount[r] += 1
    return spCount

def calc_twod_spectra(u):
    ny,nx = np.shape(u)
    umean, ufluct = get_mean_fluct(u)
    nWavenumbers = np.int(nx/np.sqrt(2)) #May not work when nx != ny
    uspectra = np.zeros(nWavenumbers)
    spCount = get_twod_sp_count(nx, ny, nWavenumbers)
    k1 = np.fft.fftfreq(nx)*nx
    k2 = np.fft.fftfreq(ny)*ny
    uhatN = np.fft.fft2(ufluct)/np.double(nx*ny)
    for i in range(nx):
        for j in range(ny):
            r = np.rint(np.sqrt(k1[i]*k1[i] + k2[j]*k2[j]))
            if(r < nWavenumbers):
                uspectra[r] += uhatN[i,j]*np.conj(uhatN[i,j])       
    
    for i in range(1,nWavenumbers):
        uspectra[i] = uspectra[i] *  np.pi * ((i+0.5)**2 - (i-0.5)**2) / spCount[i]
    
    return uspectra


def calc_L(utau, g, To, Qo, kappa):
    L = np.empty(np.shape(utau))
    L = -utau * utau * utau / (kappa * g * Qo / To)
    return L

def calc_zi(wpthetap, zvar):
    nt, nz = np.shape(wpthetap)
    zi = np.empty(nt)    
    for it in range(nt):
        zi[it] = zvar[np.argmin(wpthetap[it,:])]
        
        
