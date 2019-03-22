import os, sys
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib/')

import dd_20180130
import sftth


toth_d = {'equ_exp': 'AUGD', 'equ_dia':'GQH', 'equ_ed': 0, \
          'out_exp': 'AUGD', 'NBIpar':'TRANSP 2012', 't_fringe': 0.,\
          'ne_exp': 'AUGD', 'ne_diag': 'DCK', 'ne_sig': 'H-0', 'ne_ed': 0}

lastshot = dd_20180130.LastShotNr()

# After 19:00, once
   
for nshot in range(lastshot-25, lastshot):
    toth_d['shot'] = nshot
    sftth.write_tot_tth(toth_d)
