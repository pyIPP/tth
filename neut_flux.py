import sys
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')
import numpy as np
import dd_20140407

sf = dd_20140407.shotfile()
ioc_sig = ('F04', 'F09', 'F03', 'F01')


class IOC_FLUX:


    def __init__(self, nshot):

        nioc = len(ioc_sig)
        ntioc = 0
        jioc = 0

        if sf.Open('IOC', nshot):
            self.tioc = sf.GetTimebase('Timebase')
            ntioc = len(self.tioc)
            self.iflux = 0
            self.neut_flux = np.empty(ntioc)
            self.neut_flux[:] = np.nan
            for sig in ioc_sig:
                trace = sf.GetSignal(sig)
                if trace is not None:
                    no_nan = np.isfinite(trace)
                    self.neut_flux[no_nan] = trace[no_nan]
                    self.iflux = jioc
                    self.ioc_flag = 0
                jioc+= 1
            sf.Close()
        else:
            self.iflux = 8

        if not hasattr(self, 'neut_flux'):
            self.ioc_flag = 1
        else:
            self.neut_flux = np.nan_to_num(self.neut_flux)

        self.iflux += 1


if __name__ == '__main__':

    import matplotlib.pylab as plt

#    nshot = 18863
    nshot = 28053
 
    ioc = IOC_FLUX(nshot)
    ntioc = len(ioc.tioc)

    dioc = {}

    if sf.Open('IOC', nshot):
        for sig in ioc_sig:
            tmp = sf.GetSignal(sig)
            if tmp is not None:
                if len(tmp) == ntioc:
                    dioc[sig] = tmp
        sf.Close()

    plt.figure(1)
 
    for sig in ioc_sig:
        if sig in dioc.iterkeys():
            trace = dioc[sig]
            plt.plot(ioc.tioc, trace, label=sig)

    if hasattr(ioc, 'neut_flux'):
        plt.plot(ioc.tioc, ioc.neut_flux, label='NEUT_FLUX*1.01')
    plt.legend()
    plt.show()
