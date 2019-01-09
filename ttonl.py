import os, sys
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib/')
sys.path.append('/afs/ipp/home/g/git/python/repository')
import sfdiff, any_shot_today
import time, subprocess
import dd_20180130
import ww_20180130
import exec_toth

sf  = dd_20180130.shotfile()

firstshot = False
# On a non-shotday the loop is terminated by the crontab final time (7 p.m.)
while not firstshot:
    firstshot = any_shot_today.anyshot()
    if not firstshot:
        time.sleep(60)

toth_d = {'equ_exp': 'AUGD', 'equ_dia':'GQH', 'equ_ed': 0, \
          'out_exp': 'AUGD', 'NBIpar':'TRANSP 2012', 't_fringe': 0.,\
          'ne_exp': 'AUGD', 'ne_diag': 'DCK', 'ne_ed': 0}

if any_shot_today.anyshot():

    hour = 0
    while (hour < 19):

        loctime = time.localtime(time.time())
        hour = loctime[3]
        print('Waiting for next shot')
        print('ctrl+c to terminate the script')
        sfhdir = os.path.dirname(os.path.realpath(__file__))
        try:
            nshot = dd_20180130.LastShotNr()
        except:
            print('Problems reading last shot number')
            time.sleep(60)
            continue

        print('Last shot: %d' %nshot)
        print('Computing TOT for shot %d' %nshot)
        toth = exec_toth.ex_toth(nshot, toth_d)
        if 'time' in toth.tot.iterkeys():
            if sfdiff.sfdiff(nshot, toth_d['out_exp'], 'TOT', toth.tot):
                    ww_20180130.write_sf(nshot, toth.tot, sfhdir, 'TOT', exp=toth_d['out_exp'])

        nshot1 = nshot - 1
        print('Computing TTH for shot %d' %nshot1)
        toth_d['shot'] = nshot1
        toth1 = exec_toth.ex_toth(nshot1, toth_d)
        if 'time' in toth1.tot.iterkeys():
            if sfdiff.sfdiff(nshot1, toth_d['out_exp'], 'TOT', toth1.tot):
                ww_20180130.write_sf(nshot1, toth1.tot, sfhdir, 'TOT', exp=toth_d['out_exp'])
            toth1.tth['TOT_file']['edition'] = -1
            if sf.Open('TOT', nshot1):
                print('\nTOT edition used for TTH: %d\n' %sf.edition)
                toth1.tth['TOT_file']['edition'] = sf.edition
                sf.Close()
            if sfdiff.sfdiff(nshot1, toth_d['out_exp'], 'TTH', toth1.tth):
                ww_20180130.write_sf(nshot1, toth1.tth, sfhdir, 'TTH', exp=toth_d['out_exp'])

        time.sleep(60)

# After 19:00, once

    for shot in range(nshot-25, nshot):
        toth_d['shot'] = shot
        exec_toth.ex_toth(toth_d)

else:
    print('Today is no shotday')
