import os, sys
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib/')
sys.path.append('/afs/ipp/home/g/git/python/repository')
import any_shot_today
import time
import dd_20180130
import sftth

firstshot = False
# On a non-shotday the loop is terminated by the crontab final time (7 p.m.)
while not firstshot:
    firstshot = any_shot_today.anyshot()
    if not firstshot:
        time.sleep(60)

toth_d = {'equ_exp': 'AUGD', 'equ_dia':'GQH', 'equ_ed': 0, \
          'out_exp': 'AUGD', 'NBIpar':'TRANSP 2012', 't_fringe': 0.,\
          'ne_exp': 'AUGD', 'ne_diag': 'DCK', 'ne_sig': 'H-0', 'ne_ed': 0}

if any_shot_today.anyshot():

    hour = 0
    while (hour < 19):

        loctime = time.localtime(time.time())
        hour = loctime[3]
        print('Waiting for next shot')
        print('ctrl+c to terminate the script')
        try:
            lastshot = dd_20180130.LastShotNr()
            flag = True
        except:
            print('Problems reading last shot number')
            flag = False
            continue

        if flag:
            for nshot in range(lastshot, lastshot-3, -1):
                toth_d['shot'] = nshot
                w_tth = (nshot < lastshot)
                try:
                    sftth.write_tot_tth(toth_d, w_tth=w_tth)
                except:
                    print('Problems executing TOTH')

        time.sleep(60)

else:
    print('Today is no shotday')
