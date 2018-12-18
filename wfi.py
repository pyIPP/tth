import numpy as np

def wfi(twfi, tau_sd, dt):

    print('\nCalculating fast ion energy')
    nttot = len(twfi)
    (wfi_pos, ) = np.where(twfi > 0)
    twdelay = np.zeros(nttot, dtype=np.float32)
    xtwfi   = np.zeros(nttot, dtype=np.float32)
    dtx = 2.6*dt

    if len(wfi_pos) == 0:
        print('TOT time base has no intersection with NBI on phase')
    else:
        wfi_beg = max(np.min(wfi_pos), 1)
        wfi_end = min(np.max(wfi_pos)+200, nttot)
        for jt in range(wfi_beg, wfi_end):
            if (tau_sd[jt] > 0) and (tau_sd[jt-1] > 0):
                twdelay[jt] = twdelay[jt-1]+0.5*dtx/tau_sd[jt]+0.5*dtx/tau_sd[jt-1]
            else:
                twdelay[jt] = twdelay[jt-1]
# d/dt twdelay = - 2.6    1/tau_sd 
# twfi is a step function, xtwfi the integral with exp decay
            (ind, ) = np.where(twdelay-twdelay[jt] < 10)
            wfi1 = max(wfi_beg, np.min(ind))
            wfi2 = min(wfi_end+1, jt)
# Trapezoidal sum: extremes
            try:
                xtwfi[jt] = 0.5*twfi[wfi1]/tau_sd[wfi1]*np.exp(twdelay[wfi1] - twdelay[jt]) + \
                            0.5*twfi[wfi2]/tau_sd[wfi2]*np.exp(twdelay[wfi2] - twdelay[jt])
                wfitmp = twfi[wfi1+1: wfi2]/tau_sd[wfi1+1: wfi2] * \
                         np.exp(twdelay[wfi1+1: wfi2] - twdelay[jt])
            except:
                print('\nwfi2 = %d\n' %wfi2)

            xtwfi[jt] += np.sum(wfitmp)

    return dtx*xtwfi
