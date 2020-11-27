import numpy as np
import dd_20180130

sf = dd_20180130.shotfile()
ioc_sig = ('F04', 'F09', 'F03', 'F01')


class IOC_FLUX:


    def __init__(self, nshot, verb=False):

        nioc = len(ioc_sig)
        ntioc = 0

        if sf.Open('IOC', nshot):
            self.tioc = sf.GetTimebase('Timebase')
            ntioc = len(self.tioc)
            self.iflux = 0
            self.neut_flux = np.empty(ntioc)
            self.neut_flux[:] = np.nan
            for jioc, sig in enumerate(ioc_sig):
                trace = sf.GetSignal(sig)
                if (trace is None) or (np.sum(trace) == 0.):
                    continue
                trace = trace[:ntioc]
                no_nan = np.isfinite(trace)
                if verb:
                    print(sig, len(no_nan), np.sum(no_nan), np.sum(trace), len(trace), len(self.neut_flux), len(self.tioc))
                for j, val in enumerate(trace):
                    if no_nan[j]:
                        self.neut_flux[j] = val
                self.iflux = jioc
                self.ioc_flag = 0
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
    import sys
    sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')

#    nshot = 18863
    nshot = 33804
 
    ioc = IOC_FLUX(nshot, verb=True)
    ntioc = len(ioc.tioc)

    dioc = {}

    if sf.Open('IOC', nshot):
        for sig in ioc_sig:
            tmp = sf.GetSignal(sig)[:ntioc]
            if tmp is not None:
                if len(tmp) == ntioc:
                    dioc[sig] = tmp
            else:
                print('Signal %s not found' %sig)
        sf.Close()

    plt.figure(1)
 
    for sig in ioc_sig:
        if sig in dioc.keys():
            trace = dioc[sig]
            plt.plot(ioc.tioc, trace, label=sig)

    if hasattr(ioc, 'neut_flux'):
        plt.plot(ioc.tioc, ioc.neut_flux, label='NEUT_FLUX*1.01')
    plt.legend()
    plt.show()
