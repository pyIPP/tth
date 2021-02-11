import dd_20200525
import numpy as np

sf  = dd_20200525.shotfile()

#-----------------------
def sfdiff(nshot, exp, diag, dcheck, toler_frac=0.0001):

# Check previous shotfile edition

    print('Checking last edition of %s:%s for shot %d' %(exp, diag, nshot) )
    dcheck_old = None
    sf_diff = False
    if sf.Open(diag, nshot, experiment=exp, edition=0):
        for sig in sf.GetNames():
            info = sf.GetInfo(sig)
            if info.error != 0:
                print('Signal %s not found in previous edition' %sig)
                return True
            otype = info.objtyp
            if otype not in (7, ):
                continue
            print('Signal %s' %sig)
            dcheck_old = sf.GetSignal(sig)
            if dcheck_old is None and dcheck[sig] is not None:
                print('Problem getting signal %s' %sig)
                return True
            if dcheck[sig].shape != dcheck_old.shape:
                print('Shape in old/new edition', \
                      dcheck_old.shape, dcheck[sig].shape )
                return True
            if np.count_nonzero(np.isnan(dcheck_old)) != np.count_nonzero(np.isnan(dcheck[sig])):
                print('Different amount of Nan')
                return True
            max_val = np.max(np.abs(dcheck[sig]))
            if max_val > 0:
                toler_val = toler_frac*max_val
                tmp = np.abs(dcheck_old - dcheck[sig])
                if (np.max(tmp) > toler_val):
                    (ind, ) = np.where(tmp > toler_val)
                    print('Difference in signal %s' %sig)
                    print(dcheck_old[ind], dcheck[sig][ind])
                    return True
        sf.Close()
    else:
        print('Problems opening %s:%s for shot %d' %(exp, diag, nshot) )
        return True

    if not sf_diff:
        print('\nSuccessful check, no new %s:%s edition written\n' %(exp, diag))

    return sf_diff
