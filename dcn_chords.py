import sys
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')

import numpy as np
import map_equ_20161123

eqm = map_equ_20161123.equ_map()

HFS_R = np.array([1.0525, 1.006, 1.0067, 1.008, 1.1287, 1.0946])
HFS_z = np.array([0.2233, 0.1447, 0.3154, -0.146, 1.0566, 0.8032])
LFS_R = np.array([2.2491, 2.1664, 2.1695, 2.1614, 2.1668, 2.1715])
LFS_z = np.array([-0.0217, 0.1533, 0.3236, -0.1603, 0.1528, 0.4426])

Rmid = 0.5*(HFS_R + LFS_R)
zmid = 0.5*(HFS_z + LFS_z)
phi_up = 180./np.pi*np.arctan((LFS_z - HFS_z)/(LFS_R -HFS_R))

n_dcn = len(Rmid)

def chord_len(nshot, diag='EQH', dcn_ch=None):

    dtot = {}
    nttot = 10000
    dt = 0.001
    dtot['time'] = dt*(1 + np.arange(nttot))
    equ_diag = diag
    eqm_status = eqm.Open(nshot, diag=diag)

    if dcn_ch == None:
        dcn_channels = range(n_dcn)
    else:
        dcn_channels = np.atleast_1d(dcn_ch)

    while not eqm_status:
        for equ_diag in ('EQH', 'EQI', 'FPP'):
            eqm_status = eqm.Open(nshot, diag=equ_diag)
            if eqm_status:
                break
    if not eqm_status:
        dtot = np.zeros((nttot, n_dcn))
        return dtot

    for jdcn in dcn_channels:
        print('Determining chord length of H-%d using equ %s' %(jdcn, diag) )
        R1, z1 = eqm.cross_surf(0.9999, r_in=Rmid[jdcn], z_in=zmid[jdcn], theta_in=np.radians(phi_up[jdcn]), coord_in='rho_pol')

        chord = np.hypot(R1[:, 0, 1] - R1[:, 0, 0], z1[:, 0, 1] - z1[:, 0, 0])

        dtot[jdcn] = np.interp(dtot['time'], eqm.t_eq, chord)
    dtot['diag'] = equ_diag

    return dtot


if __name__ == '__main__':

    import matplotlib.pylab as plt
    import dd_20140407

    sf = dd_20140407.shotfile()
#    dtot = cut_surf(14481, diag='FPP')

    shot = 19430
#    shot = 28053
    if False:
        if eqm.Open(shot, diag='EQH'):
            eqm.read_ssq()
            for sig in ('H-1', 'H-2', 'H-3', 'H-4', 'H-5', 'V-1', 'V-2'):
                plt.plot(eqm.t_eq, eqm.ssq['len%s' %sig], label=sig)
            plt.show()

    if True:
        dtot = cut_surf(shot, diag='EQH')
        for key in R.keys():
            plt.plot(dtot['time'], dtot['len%s' %key], label=key)
        plt.figtext(0.5, 0.85, dtot['diag'], ha='center')
        plt.legend()
        plt.show()
