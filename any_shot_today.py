import sys, datetime
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')

import dd_20180130

sf = dd_20180130.shotfile()


def anyshot():

    nshot = dd_20180130.LastShotNr()

    today = datetime.datetime.now().strftime('%d%b%Y')
    date_last = 'No date'
    if sf.Open('JOU', nshot):
        date_last = sf.date[:9].strip()
        sf.Close()
    try:
        day = int(date_last[:2])
    except:
        date_last = '0%s' %date_last
    print('Today is     %s' %today) 
    print('Last shot on %s' %date_last) 

    return(today == date_last)

if __name__ == '__main__':

    print(anyshot())
