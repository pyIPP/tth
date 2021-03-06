import dd_20180130 as dd

sf = dd.shotfile()

def  extract_main(gasv):

    gasv_str = gasv.decode('utf-8') if isinstance(gasv, bytes) else gasv
    main_spec = 'D'
    if '(1)' in gasv_str:
        main_spec = gasv_str.split('(1)')[0].split(':')[-2].split(',')[-1].strip()
    elif gasv_str.count(':') == 1:
        main_spec = gasv_str.split(':')[0].strip()
    return main_spec

def spec_jou_sf(nshot):

    if sf.Open('JOU', nshot):
        tmp = sf.GetParameter('FILLING', 'GasVent')
        sf.Close()
        gasv_str = b''.join(tmp)
        return extract_main(gasv_str)
    else:
        return None
