import sys
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib/')
import numpy as np
from sf2equ_20200525 import READ_EQU
import mapeq_20200507 as meq


HFS_R = np.array([1.0525, 1.006, 1.0067, 1.008, 1.1287, 1.0946])
HFS_z = np.array([0.2233, 0.1447, 0.3154, -0.146, 1.0566, 0.8032])
LFS_R = np.array([2.2491, 2.1664, 2.1695, 2.1614, 2.1668, 2.1715])
LFS_z = np.array([-0.0217, 0.1533, 0.3236, -0.1603, 0.1528, 0.4426])

Rmid = 0.5*(HFS_R + LFS_R)
zmid = 0.5*(HFS_z + LFS_z)
phi_up = np.degrees(np.arctan2(LFS_z - HFS_z, LFS_R - HFS_R))

n_dcn = len(Rmid)

def chord_len(nshot, diag='EQH', dcn_ch=None):

    dtot = {}
    nttot = 10000
    dt = 0.001
    dtot['time'] = dt*(1. + np.arange(nttot, dtype=np.float32))

    if dcn_ch == None:
        dcn_channels = range(n_dcn)
    else:
        dcn_channels = np.atleast_1d(dcn_ch)

    eqsf = READ_EQU(nshot, diag='EQH', exp='AUGD', ed=0)
    if not hasattr(eqsf, 'shot'):
        eqsf = READ_EQU(nshot, diag='EQI', exp='AUGD', ed=0)
        if not hasattr(eqsf, 'shot'):
            eqsf = READ_EQU(nshot, diag='FPP', exp='AUGD', ed=0)
            if not hasattr(eqsf, 'shot'):
                return np.zeros( (nttot, len(dcn_channels)) )
    eqsf.read_scalars()
    eqsf.read_profiles()
    eqsf.read_pfm()
    eqsf.close()

    for jdcn in dcn_channels:
        if (nshot < 32000) and (jdcn == 0):
            dtot[jdcn] = np.ones(nttot, dtype=np.float32)
        else:
            print('Determining chord length of H-%d using equ %s' %(jdcn, diag) )
            R1, z1 = meq.cross_sep(eqsf, r_in=Rmid[jdcn], z_in=zmid[jdcn], theta_in=np.radians(phi_up[jdcn]) )

            chord = np.hypot(R1[:, 1] - R1[:, 0], z1[:, 1] - z1[:, 0])

            dtot[jdcn] = np.interp(dtot['time'], eqsf.time, chord) + 1.e-8
            print('First intersection R: %8.4f' %R1[0, 0] )
            del R1, z1
    dtot['diag'] = eqsf.diag
    dtot['ed']   = eqsf.ed

    return dtot


if __name__ == '__main__':

    import matplotlib.pylab as plt
    import dd_20200525 as dd

    sf = dd.shotfile()

    shot = 34027

    dtot = chord_len(shot, diag='EQH', dcn_ch=0)
    plt.plot(dtot['time'], dtot[0])
    plt.figtext(0.5, 0.85, dtot[0], ha='center')
    plt.show()
