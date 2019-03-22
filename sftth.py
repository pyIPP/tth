import os, sys
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib/')
sys.path.append('/afs/ipp/home/g/git/python/repository')
import sfdiff
import dd_20180130, ww_20180130
import exec_toth

sf = dd_20180130.shotfile()

def write_tot_tth(toth_d, w_tth=True):

    nshot = toth_d['shot']
    try:
        toth = exec_toth.ex_toth(nshot, toth_d)
    except:
        print('TOT/TTH not executed')
        return
    sfhdir = os.path.dirname(os.path.realpath(__file__))
    if 'time' in toth.tot.iterkeys():
        if sfdiff.sfdiff(nshot, toth_d['out_exp'], 'TOT', toth.tot):
            ww_20180130.write_sf(nshot, toth.tot, sfhdir, 'TOT', exp=toth_d['out_exp'])

        toth.tth['TOT_file']['expr']    = toth_d['out_exp']
        if sf.Open('TOT', nshot, experiment=toth_d['out_exp']):
            print('TOT edition used for TTH: %d' %sf.edition)
            toth.tth['TOT_file']['edition'] = sf.edition
            sf.Close()

            if w_tth and sfdiff.sfdiff(nshot, toth_d['out_exp'], 'TTH', toth.tth):
                ww_20180130.write_sf(nshot, toth.tth, sfhdir, 'TTH', exp=toth_d['out_exp'])
