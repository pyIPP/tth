import sys, struct
import numpy as np
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')
from sfarr import SFARR


def bin2sfarr(f, shape=None, fmt='f'):

    struct_len = struct.calcsize(fmt)
    if shape is None:
        return struct.unpack(fmt, f.read(struct_len))[0]
    nlen = np.prod(shape)
    if fmt == 's':
        return str(struct.unpack('%ds' %nlen, f.read(nlen))[0])
    tmp = struct.unpack('%d%s' %(nlen, fmt), f.read(struct_len*nlen))
    arr = np.array(tmp).reshape(shape)

    return SFARR(np.nan_to_num(arr))
