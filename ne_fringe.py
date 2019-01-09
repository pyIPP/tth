import numpy as np
import dd_20140407

sf  = dd_20140407.shotfile()

dcn_sig = { 'DCS': 'nedl_H1', 'DCK': 'H-1', \
            'DCN': 'H-1', 'DCR': 'H-1', 'DCP': 'H-1'}

def fringe_jumps(time_in, sig_in, flattop_end):

# Fringe jump detection

# delta_t = 0.25ms <-> interferometer signal filter at 4 kHz
# delta_n = 0.572e19 <-> 1 fringe

    if type(time_in) == type(0.1):
        print('Time base is a single point, skipping fringe jump detection')
        return 0.0
    delta_n_fj = 0.5e19
    delta_t_filter = 0.25e-3
    print('Max time_in = %8.4f,  flattop = %8.4f' %(np.max(time_in), flattop_end) )
    if (np.max(time_in) <= flattop_end):
        print('No time points after flat top end, skipping fringe jump detection')
        return flattop_end
    jtflattop_end = np.min(np.where(time_in > flattop_end))        
    jtdiff = int(delta_t_filter/(time_in[1]-time_in[0])) + 1
    jtmin = np.min(np.where(time_in > 0.5))
# git    jtmin = max(5, jtdiff+1)
    jtjump = len(time_in)-1
    for jt in range(jtmin, jtflattop_end):
        sig_t1 = sig_in[jt-jtdiff]
        sig_t2 = sig_in[jt]
        if (sig_t2 < 0):
            print('Negative density value')
            return time_in[jt]
        else:
            if (abs(sig_t1-sig_t2) > delta_n_fj):
                jtjump = jt
                break
    return time_in[jtjump]


def dc_jump(nshot, diag_in, flattop_end=10., exp_in='AUGD', ed_in=0):

    tjump_dc = 0
    ne_ed = None

    if diag_in not in dcn_sig.keys():
        print('Diag %s not supported' %diag_in)
        return tjump_dc, ne_ed

    if sf.Open(diag_in, nshot, experiment=exp_in, edition=ed_in):
        ne_sig = dcn_sig[diag_in]
        tdcn = sf.GetTimebase(ne_sig)
        dcn_ref = sf.GetSignal(ne_sig, cal=True)
        ne_ed = sf.edition
        sf.Close()

        tjump_dc = fringe_jumps(tdcn, dcn_ref, flattop_end)
        print('Fringe jump at %8.4f' %tjump_dc )
        print('Shot %s, %s:%s' %(nshot, diag_in, ne_sig) )

    else:
        print('Shotfile %s not found for #%d' %(diag_in, nshot))

    return tjump_dc, ne_ed


def ne_fringe(nshot, flattop_end=10., exp_in='AUGD', diag_in='DCK', ed_in=0, tj_forced=None):

# Hierarchy: 
#1) prescribed (if any)
#2) DCK
#3) DCN

    print('\nDensity signal')

    ne_diag = diag_in
    ne_exp  = exp_in
    ne_ed   = ed_in
    tjump_dc, dcn_ed = dc_jump(nshot, ne_diag, flattop_end=flattop_end, exp_in=ne_exp, ed_in=ne_ed)
    if (dcn_ed is not None) and (tj_forced is not None):
        tjump_dc = tj_forced

#AUGD:DCK(0) 1st fallback

    if dcn_ed is None:
        ne_exp  = 'AUGD'
        ne_diag = 'DCK'
        ne_ed   = 0
        tjump_dc, dcn_ed = dc_jump(nshot, ne_diag, flattop_end=flattop_end, exp_in=ne_exp, ed_in=ne_ed)

#AUGD:DCN(0) 2nd fallback

    if dcn_ed is None:
        ne_exp  = 'AUGD'
        ne_diag = 'DCN'
        ne_ed   = 0
        tjump_dc, dcn_ed = dc_jump(nshot, ne_diag, flattop_end=flattop_end, exp_in=ne_exp, ed_in=ne_ed)

# What was really used

    if dcn_ed is None: # no DCN found
        fringe_d = None
    else:
        tjump = tjump_dc
        fringe_d = {'exp': ne_exp, 'diag': ne_diag, \
                    'sig': dcn_sig[ne_diag], 'ed': dcn_ed, 'tjump': tjump}

    return fringe_d
