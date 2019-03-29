import numpy as np
import map_equ_20161123

eqm = map_equ_20161123.equ_map()

#Rcoo = [1.785, 1.20]
#R1dcn = [1.008, 1.007, 1.008, 1.12869, 1.09459]
#z1dcn = [0.146, 0.317, -0.146, 1.05664, 0.803226]
#R2dcn = [2.1664, 2.1695, 2.1614, 2.1668, 2.1715]
#z2dcn = [0.1533, 0.3236, -0.1603, 0.1528, 0.4426]

#Rmid = 0.5*(R1dcn + R2dcn)
#zmid = 0.5*(z1dcn + z2dcn)
#phi_up = 180.0/np.pi*np.arctan((z2dcn[jdcn]-z1dcn[jdcn])/(R2dcn[jdcn]-R1dcn[jdcn]))

Rmid = {'V-1': 1.7850, 'V-2': 1.20, \
     'H-1': 1.5872, 'H-2': 1.5882, 'H-3': 1.5847  , 'H-4': 1.6477, 'H-5': 1.633}
zmid = {'V-1': 0., 'V-2': 0., \
     'H-1':0.14965, 'H-2': 0.3203, 'H-3': -0.15315, 'H-4': 0.60472, 'H-5': 0.62291}
phi_up = {'V-1': 90., 'V-2': 90., 'H-1': 0.36106, \
          'H-2': 0.3253, 'H-3': -0.71032, 'H-4': -41.045, 'H-5': -18.514}


def cut_surf(nshot, exp='AUGD', diag='EQH', sig=None, min_len=0., max_len=2.5):

    dtot = {}
    nttot = 10000
    dt = 0.001
    dtot['time'] = dt*(1 + np.arange(nttot))
    equ_diag = diag
    eqm_status = eqm.Open(nshot, diag=diag)

    if sig is None:
        R_list = Rmid.keys()
    else:
        R_list = [sig]

    while not eqm_status:
        for equ_diag in ('EQH', 'EQI', 'FPP'):
            eqm_status = eqm.Open(nshot, diag=equ_diag)
            if eqm_status:
                break
    if not eqm_status:
        for sig in R_list:
            dtot[lensig] = np.zeros(nttot)
        return dtot

    for sig in R_list:
        lensig = 'len%s' %sig
        R1, z1 = eqm.cross_surf(0.999, r_in=Rmid[sig], z_in=zmid[sig], theta_in=np.radians(phi_up[sig]), coord_in='rho_pol')

        chord = np.hypot(R1[:, 0, 1] - R1[:, 0, 0], z1[:, 0, 1] - z1[:, 0, 0])

        chord = np.maximum(chord, min_len)
        chord = np.minimum(chord, max_len)

        dtot[lensig] = np.interp(dtot['time'], eqm.t_eq, chord)
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
