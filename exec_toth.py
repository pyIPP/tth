import os
import numpy as np
import scipy.ndimage
import NBIlosses, tderiv, wfi, ne_fringe, dcn_chords, ne_sep_sol
import dd_20140407

sf  = dd_20140407.shotfile()


parlaws = { \
    'prefac':'coeff', 'a_i':'A', 'Ip':'IpiFP', 'P_tot':'P_TOT', 'P_net':'P_NET', \
    'Rgeo':'Rgeo', 'ahor':'ahor', 'kappa':'kappa', 'karea':'karea', \
    '<ne>_H-1':'H-1_corr', 'Btf':'BTF' \
                }

laws = { \
    'ITERL-89P(tot), Wtot/Ptot': \
        {'scalmod': 'tau_tot', 'coeff':0.038, 'A':0.5, 'IpiFP':0.85, 'P_TOT':-0.5, 'P_NET':0,            \
         'Rgeo':1.2, 'ahor':0.3, 'kappa':0.5, 'karea':0, 'H-1_corr':0.1, 'BTF':0.2},            \
    'ITERL-96P(th), Wth/Pnet': \
        {'scalmod': 'tau_th', 'coeff':0.023, 'A':0.2, 'IpiFP':0.96, 'P_TOT':0, 'P_NET':-0.73,         \
         'Rgeo':1.89, 'ahor':-0.06, 'kappa':0.64, 'karea':0, 'H-1_corr':0.4, 'BTF':0.03}, \
    'ITERH92PY:tau_E, ELMy, Wtot/Pnet': \
        {'scalmod': 'tau_hybr', 'coeff':0.051, 'A':0.51, 'IpiFP':0.83, 'P_TOT':0, 'P_NET':-0.51,        \
         'Rgeo':1.76, 'ahor':0.11, 'kappa':0.5, 'karea':0, 'H-1_corr':0.05, 'BTF':0.1},     \
    'ITERH-92P(y, tot), ELMy, Wtot/Ptot': \
        {'scalmod': 'tau_tot', 'coeff':0.038, 'A':0.38, 'IpiFP':0.76, 'P_TOT':-0.5, 'P_NET':0,         \
         'Rgeo':1.87, 'ahor':0.2, 'kappa':0.8, 'karea':0, 'H-1_corr':0.11, 'BTF':0.19},     \
    'ITERH-92P(y, th), ELMy, Wth/Pnet': \
        {'scalmod': 'tau_th', 'coeff':0.034, 'A':0.4, 'IpiFP':0.90, 'P_TOT':0, 'P_NET':-0.65,         \
         'Rgeo':1.9, 'ahor':0.2, 'kappa':0.8, 'karea':0, 'H-1_corr':0.3, 'BTF':0.05},         \
    'ITERH-93P(th), ELM free, Wth/Pnet': \
        {'scalmod': 'tau_th', 'coeff':0.036, 'A':0.41, 'IpiFP':1.06, 'P_TOT':0, 'P_NET':-0.67,        \
         'Rgeo':1.9, 'ahor':-0.11, 'kappa':0.66, 'karea':0, 'H-1_corr':0.17, 'BTF':0.32}, \
    'ITERH-98P(y, th), ELMy, Wth/Pnet': \
        {'scalmod': 'tau_th', 'coeff':0.0365, 'A':0.2, 'IpiFP':0.97, 'P_TOT':0, 'P_NET':-0.63,        \
         'Rgeo':1.70, 'ahor':0.23, 'kappa':0.67, 'karea':0, 'H-1_corr':0.41, 'BTF':0.08}, \
    'ITERH-98P(y, th, 2), ELMy, Wth/Pnet': \
        {'scalmod': 'tau_th', 'coeff':0.0562, 'A':0.19, 'IpiFP':0.93, 'P_TOT':0, 'P_NET':-0.69,     \
         'Rgeo':1.39, 'ahor':0.58, 'kappa':0, 'karea':0.78, 'H-1_corr':0.41, 'BTF':0.15}, \
    'DS03 or ESGB, McDonald, ppcf 46 (2004) A215 equation 11': \
        {'scalmod': 'tau_th', 'coeff':0.028, 'A':0.14, 'IpiFP':0.83, 'P_TOT':0, 'P_NET':-0.55,        \
         'Rgeo':1.81, 'ahor':0.3, 'kappa':0.75, 'karea':0, 'H-1_corr':0.49, 'BTF':0.07},    \
    'CORDEY05, NF 45 (2005) 1078, equation 9': \
        {'scalmod': 'tau_th', 'coeff':0.0506, 'A':0.11, 'IpiFP':0.85, 'P_TOT':0, 'P_NET':-0.45,     \
         'Rgeo':1.21, 'ahor':0.39, 'kappa':0, 'karea':0.82, 'H-1_corr':0.26, 'BTF':0.17},    \
    'KARDAUN, IAEA 2006, IT/P10': \
        {'scalmod': 'tau_tot', 'coeff':0, 'A':0, 'IpiFP':0, 'P_TOT':0, 'P_NET':0,     \
         'Rgeo':0, 'ahor':0, 'kappa':0, 'karea':0, 'H-1_corr':0, 'BTF':0} \
#    'KARDAUN-LANG, IAEA 2012, EX/P4-01': \
#        {'scalmod': 'tau_tot', 'coeff':0, 'A':0, 'IpiFP':0, 'P_TOT':0, 'P_NET':0,     \
#         'Rgeo':0, 'ahor':0, 'kappa':0, 'karea':0, 'H-1_corr':0, 'BTF':0}
}

tth_laws = ( \
    'ITERL-89P(tot), Wtot/Ptot', \
    'ITERL-96P(th), Wth/Pnet', \
    'ITERH92PY:tau_E, ELMy, Wtot/Pnet', \
    'ITERH-92P(y, tot), ELMy, Wtot/Ptot', \
    'ITERH-92P(y, th), ELMy, Wth/Pnet', \
    'ITERH-93P(th), ELM free, Wth/Pnet', \
    'ITERH-98P(y, th), ELMy, Wth/Pnet', \
    'ITERH-98P(y, th, 2), ELMy, Wth/Pnet', \
    'DS03 or ESGB, McDonald, ppcf 46 (2004) A215 equation 11', \
    'CORDEY05, NF 45 (2005) 1078, equation 9', \
    'KARDAUN, IAEA 2006, IT/P10' )
#     'KARDAUN-LANG, IAEA 2012, EX/P4-01' \


#-----------------------
def msg_quit(nshot, msg):

    import datetime

    msg2 = str(nshot) + ' '+msg+'\n'
    logfile = 'toth.log'

    inlog = open(logfile, 'r')
    lines = inlog.readlines()
    inlog.close()

    if msg2 not in lines:
        flog = open(logfile, 'a')
        print('Program stopped')
        print(msg)
        now = datetime.datetime.now()
        tout = now.strftime("%d-%m-%Y %H:%M:%S\n")
        flog.write(tout)
        flog.write(msg2)
        flog.close()
    return


#-----------------------
def sig2toth(ttot, nshot, diag, sig, exp='AUGD', ed=0, min_zero=0.01):

# Read signal and interpolate to TOT time base

    print('Reading %s from diag %s and mapping it onto TOT time base' %(sig, diag) )

    d_out = None
    if sf.Open(diag, nshot, experiment=exp, edition=ed):
        info = sf.GetInfo(sig)
        if info.error == 0:
            tdiag = sf.GetTimebase(sig)
            if tdiag is None:
                return d_out
            diag_sig = sf.GetSignal(sig, cal=True)
            print np.max(diag_sig), len(diag_sig), len(tdiag)
            if len(tdiag) == len(diag_sig):
                diag_sig = np.nan_to_num(diag_sig)
                d_out    = np.interp(ttot, tdiag, diag_sig, left=0, right=0)
                ind_zero = np.where(d_out == 0)
                d_out[ind_zero] = min_zero
        sf.Close()

    return d_out


#-------------
class ex_toth:
#-------------


    def __init__(self, nshot, toth_d):


        print('\nStarting TOT/TTH evaluation')
        try:

            equexp  = toth_d['equ_exp'].strip()
            equdiag = toth_d['equ_dia'].strip()
            equed   = int(toth_d['equ_ed'])
            neexp   = toth_d['ne_exp'].strip()
            nediag  = toth_d['ne_diag'].strip()
            need    = int(toth_d['ne_ed'])
            exp_write = toth_d['out_exp'].strip()
            NBIpar = toth_d['NBIpar'].strip()[:6]
            ion_mass = 2
            t_fr = float(toth_d['t_fringe'])
            if t_fr == 0.:
                t_fr = None

            if 'spec' in toth_d.keys():
                if toth_d['spec'].strip() == 'H':
                    ion_mass = 1
                if toth_d['spec'].strip() == 'He':
                    ion_mass = 4

            home_dir = os.getenv('HOME')
            os.system('mkdir -p '+home_dir+'/shotfiles/TOT')
            os.system('mkdir -p '+home_dir+'/shotfiles/TTH')

            n_source = {}
            n_source['NBI'] = 8
            n_source['ICRH'] = 4
            n_source['ECRH'] = 8

# TOT parameter sets

            totpar = ('heat_par', 'Setup', 'Equ', 'Fringe')

# TTH parameter sets

            tthpar = ('loss_par', 'TOT_file', 'scal_par', 'G_gw_par', 'L2H_par')

            self.tot = {}
            self.tth = {}
            for pn in totpar:
                self.tot[pn] = {}
            for pn in tthpar:
                self.tth[pn] = {}

# Initial time grid
            nt_init = 10000
            dt = 0.001
            ttot_full = dt*(1 + np.arange(nt_init, dtype=np.float32))

# ====================== 
# Fetch equilibrium data
# ====================== 

#EQU


            if sf.Open(equdiag, nshot, experiment=equexp, edition=equed):
                print( 'Equilibrium data for shot %d %s:%s(%d)' \
                        %(nshot, equexp, equdiag, sf.edition) )
            else:
                if equdiag == 'FPG':
                    print('%s:FPG(%d) not found, using AUGD:GQH(0) instead' %(equexp, equed))
                    equdiag = 'GQH'
                elif equdiag == 'GQH':
                    print('%s:GQH(%d) not found, using AUGD:FPG(0) instead' %(equexp, equed))
                    equdiag = 'FPG'
# Fallback
                if sf.Open(equdiag, nshot, experiment='AUGD', edition=0):
                    print( 'Equilibrium data for shot %d AUGD:%s(%d)' \
                        %(nshot, equdiag, sf.edition) )
                    pass
                else:
                    msg = 'Neither GQH nor FPG available'
                    msg_quit(nshot, msg)
                    return
            self.tot['Equ']['exp']  = equexp
            self.tot['Equ']['diag'] = equdiag
            self.tot['Equ']['ed']   = sf.edition

            print('Equilibrium data from diag %s' %equdiag)
            tequ = sf.GetTimebase('TIMEF')
            dequ = {}

            for sig in ('Wmhd', 'Vol', 'q95', 'Rgeo', 'Raus', 'Rin', 'Circumf', 'Zoben', 'Zunt', 'fbnd-f12'):
                dequ[sig] = sf.GetSignal(sig)
            sf.Close()
            rgeo = 0.5*(dequ['Raus'] + dequ['Rin'])
            ahor = 0.5*(dequ['Raus'] - dequ['Rin'])

            den_min = 0.01
            a1 = np.maximum(ahor, den_min)
            R1 = np.maximum(rgeo, den_min)
            area = dequ['Vol']/(2.*np.pi*R1)
            kappa = 0.5*(dequ['Zoben'] - dequ['Zunt'])/a1
            karea = area/(np.pi*a1**2)
            sur = 2*np.pi*rgeo*dequ['Circumf']

#-------------------
# Flattop evaluation
#-------------------

            print('\nBtf, Ipl for flattop end estimate')

            diag = 'FPC'
            if sf.Open(diag, nshot):
                sig  = 'IpiFP'
                ipl = sf.GetSignal(sig)
                tfpc = sf.GetTimebase(sig)

            ipl_ref = 0.2*max(ipl)
            ind_down =  (ipl < ipl_ref ) & (tfpc > 1)
            flattop_end = tfpc[ind_down][-1] - 0.8
            print('Estimated flat top end: %8.4f' %flattop_end)

#-----------------------
# Fringe jumps detection
#-----------------------

            flattop_equ_end = min(tequ[-1], flattop_end)
            self.tot['Fringe'] = ne_fringe.ne_fringe(nshot, flattop_end=flattop_equ_end, \
                exp_in=neexp, diag_in=nediag, ed_in=need, tj_forced=t_fr)

            print('Time of fringe jump: %8.4f' %self.tot['Fringe']['tjump'])

#------------------
# Set TOT time base
#------------------

            t_end = min(ttot_full[-1], tequ[-1], self.tot['Fringe']['tjump'])
            print('\nT_end is min of %8.4f, %8.4f, %8.4f' %(ttot_full[-1], tequ[-1], self.tot['Fringe']['tjump']) )
            nttot = np.where(ttot_full < t_end)[0][-1]

            ttot = ttot_full[:nttot]
            self.tot['time'] = ttot
            self.tth['time'] = ttot

#----------------------
# Line averaged density
#----------------------

# Used in NBI losses, slowing-down time, H scaling laws, Schweinzer SOL scaling

            print('\nSubstract offset, normalise by chord length')

            exp  = self.tot['Fringe']['exp']
            diag = self.tot['Fringe']['diag']
            sig  = self.tot['Fringe']['sig']
            ed   = self.tot['Fringe']['ed']
            print('DC signal: %s:%s:%s(%d)' %(exp, diag, sig, ed) )

            ne_h1 = 1e-19*sig2toth(ttot, nshot, diag, sig)
            offset_h1 = np.average(ne_h1[0: 10]) # First 10 ms
            min_chord_len = 0.05
            if diag == 'DCS':
                ne_lav = ne_h1 - offset_h1 # unit 10**19 m**-3
                self.tot['H-1_corr'] = ne_lav
                self.tot['H-4_corr'] = np.zeros_like(ttot)
                self.tot['peak'] = np.ones_like(ttot)
            else:
                ne_h4 = 1e-19*sig2toth(ttot, nshot, diag, 'H-4')
                offset_h4 = np.average(ne_h4[0: 10]) # First 10 ms
                len_h1 = sig2toth(ttot, nshot, equdiag, 'lenH-1')
                if len_h1 is None:
                    len_d14 = dcn_chords.chord_len(nshot, diag='EQH', sig=['H-1', 'H-4'])
                    len_h1 = len_d14['H-1'][:nttot]
                    len_h4 = len_d14['H-4'][:nttot]
                else:
                    len_h4 = sig2toth(ttot, nshot, equdiag, 'lenH-4')
                len_h1 = np.maximum(len_h1, min_chord_len)
                len_h4 = np.maximum(len_h4, min_chord_len)
                ne_lav    = (ne_h1 - offset_h1)/len_h1
                ne_lav_h4 = (ne_h4 - offset_h4)/len_h4

                pf = ne_lav/ne_lav_h4
                pf = np.maximum(pf, 0.8)

                self.tot['H-1_corr'] = ne_lav
                self.tot['H-4_corr'] = ne_lav_h4
                self.tot['peak'] = np.minimum(pf, 4)


#-----
# Zeff
#-----

            print('Zeff')

            max_zef = 4.5
            ind_one = (ne_lav > 1)
            self.tot['Zeff'] = np.minimum(max_zef, 11.2/(ne_lav - 0.74))
            self.tot['Zeff'][~ind_one] = max_zef

#==================
# Fetch diagnostics
#==================

#----
# FPC
#----
            diag = 'FPC'
            ipl_sign = 'straight'

            for sig in ('IpiFP', 'BTF'):
                self.tot[sig] = sig2toth(ttot, nshot, diag, sig)
                if np.sum(self.tot[sig]) < 0:
                    self.tot[sig] *= -1
                    if sig == 'IpiFP':
                        ipl_sign = 'rev'
            self.tot['Setup']['ipl_sign'] = ipl_sign
#----
# JOU
#----
            jouflag = {}
            diag = 'JOU'
            if sf.Open(diag, nshot):
                for heat_lbl in ('NBI', 'ICRH', 'ECRH'):
                    jouflag[heat_lbl] = False
                    for jsrc in range(n_source[heat_lbl]):
                        pn = '%s%dL' %(heat_lbl, jsrc+1)
                        pval = sf.GetParameter(heat_lbl, pn)
                        if pval is not None:
                            if pval > 0:
                                jouflag[heat_lbl] = True
                                break
#                    pval = sf.GetParameter('FILLING', 'GAS')
#                    if not np.isnan(pval):
#                        if pval > 0:
#                            Amass = pval
                sf.Close()

            if (ion_mass != 2):
                if jouflag['NBI']:
                    msg = 'No TTH for non deuterium plasmas, NBI present'
                    msg_quit(nshot, msg)

            self.tot['Setup']['A'] = int(ion_mass)

#---
#MAG
#---

            sig = 'ULid12'
            if (nshot > 27404):
                diag = 'MAY'
                if sf.Open(diag, nshot):
                    lid12 = sf.GetSignal(sig, cal = True) 
                    tmag = sf.GetTimebase(sig) 
                    l12smo = scipy.ndimage.gaussian_filter1d(lid12, 10, axis = -1, order = 0, output = None, mode = 'reflect', cval = 0.)
                    ULid12 = np.interp(ttot, tmag, l12smo)
                    sf.Close()
            else:
                diag = 'MAG'
                if sf.Open(diag, nshot):
                    lid12 = sf.GetSignal(sig, cal = True) 
                    tmag = sf.GetTimebase(sig) 
                    ULid12 = np.interp(ttot, tmag, lid12)
                    sf.Close()
            ind_base = (tmag >= -10) & (tmag <= -9)
            if np.sum(ind_base) > 1:
                baseline = np.average(lid12[ind_base])
            else:
                baseline = np.average(lid12[:100])
            ULid12 -= baseline

            if max(np.abs(ULid12)) < 0.1:
                msg = 'No sensible Uloop data'
                msg_quit(nshot, msg)
                return
            rc = 0.05 #R = 5.e3 Ohm, C = 1.e-5 F
            xflux = np.interp(ttot, tequ, dequ['fbnd-f12'])
            self.tot['du_loop'] = tderiv.corr_u(ttot, xflux, rc)
# Used in P_OH, Schweinzer scaling
            self.tot['u_loop'] = ULid12 - self.tot['du_loop']

# ======================== 
# Map onto TOT time grid
# ======================== 

            for sig in ('Wmhd', 'q95', 'Vol'):
                self.tot[sig] = np.interp(ttot, tequ, dequ[sig])
            self.tot['Rgeo']    = np.interp(ttot, tequ, rgeo)
            self.tot['ahor']    = np.interp(ttot, tequ, ahor)
            self.tot['kappa']   = np.interp(ttot, tequ, kappa)
            self.tot['karea']   = np.interp(ttot, tequ, karea)
            self.tot['Area']    = np.interp(ttot, tequ, area)
            self.tot['Surface'] = np.interp(ttot, tequ, sur)
            self.tot['q95']     = abs(self.tot['q95'])

# ======= 
# Power
# ======= 

# Factors

            self.tth['loss_par']['fac_sol']  = 0.1
            self.tth['loss_par']['fac_orb']  = 0
            self.tot['heat_par']['fac_ecrh'] = 1

# Ohmic

            self.tot['P_OH'] = np.abs(self.tot['u_loop']*self.tot['IpiFP'])

# ICRH

            icp_min = 1e4
            if (nshot > 8218):
                sig = 'PICRFc'
                self.tot['heat_par']['fac_icrh'] = 1.0
            else:
                sig = 'PICRH'
                self.tot['heat_par']['fac_icrh'] = 0.8
            diag = 'ICP'
            picrf = sig2toth(ttot, nshot, diag, sig)
            self.tot['PICR_TOT'] = np.zeros_like(ttot)
            if picrf is not None:
                ind = (picrf >= icp_min)
                self.tot['PICR_TOT'][ind] = picrf[ind]
            if ( (nshot == 16178) or (nshot == 16201) ):
                self.tot['PICR_TOT'] *= 1.333

# ECRH

            diag = 'ECS'
            sig = 'PECRH'
            ecs_min = 1e4
            nshot_sav = nshot
            if (nshot == 15524):
                nshot = 15536
            self.tot['PECR_TOT'] = np.zeros_like(ttot)
            pecs = sig2toth(ttot, nshot, diag, sig)
            if pecs is not None:
                ind = (pecs >= ecs_min)
                self.tot['PECR_TOT'][ind] = pecs[ind]
            if (nshot_sav == 15524):
                nshot = 15524

# NBI

            diag = 'NIS'
            dnis = {}
            pow_frac = {}
            n_nbi = {}
            nnbi = 0
            min_nbi = 1.e5
            if sf.Open(diag, nshot):
                for pn in ('INJ1', 'INJ2'):
                    dnis[pn] = {}
                    for ps in ('UEXQ', 'SPEC', 'M'):
                        dnis[pn][ps] = sf.GetParameter(pn, ps)
                    n_nbi[pn] = len(dnis[pn]['UEXQ'])
                    nnbi+= len(dnis[pn]['UEXQ'])
                    n_mix = len(dnis[pn]['SPEC'])

                pniq = np.zeros((nttot, nnbi), dtype=np.float32)
                sig = 'PNIQ'
                tnis = sf.GetTimebase(sig)
                datnis = sf.GetSignalGroup(sig)
                sf.Close()
                ni_eny = np.zeros(nnbi, dtype=np.float32)
                jsrc = 0
                for jbox, pn in enumerate(('INJ1', 'INJ2')):
                    den = np.sum(dnis[pn]['SPEC'])
                    if den > 0:
                        pow_frac[pn] = dnis[pn]['SPEC'][::-1]/den
                    else:
                        pow_frac[pn] = np.zeros(n_mix, dtype=np.float32)
                    spec = (dnis[pn]['SPEC'][0]/3 + dnis[pn]['SPEC'][1]/2 + dnis[pn]['SPEC'][2])/83.2
                    print('%s power fraction: %7.4f %7.4f %7.4f' %(pn, pow_frac[pn][0], pow_frac[pn][1], pow_frac[pn][2]) )
                    for jnbi in range (n_nbi[pn]):
                        ni_eny[jsrc] = dnis[pn]['UEXQ'][jnbi]*spec
                        pniq[:, jsrc] = np.interp(ttot, tnis, datnis[:, jnbi, jbox], left = 0, right = 0)
                        jsrc += 1
                self.tot['PNBI_TOT'] = np.sum(pniq, axis=1)
            else:
                print('No NIS data found')
                self.tot['PNBI_TOT'] = np.zeros_like(ttot)

# ========== 
# dWmhd/dt
# ========== 

            print('Time derivative')
            ptst = self.tot['heat_par']['fac_icrh']*self.tot['PICR_TOT'] + \
                   self.tot['heat_par']['fac_ecrh']*self.tot['PECR_TOT'] + self.tot['PNBI_TOT']
            print('Max effective power (used to detect power jumps) %10.4e' %np.max(ptst))
            if np.max(ptst) > 0:
                self.tot['dWmhd/dt'] = tderiv.dvdt_ic(self.tot['Wmhd'], ptst)
            else:
                self.tot['dWmhd/dt'] = np.zeros_like(ttot)

# ================================== 
# Estimate ne at SOL, sep
# only reason for TOT-TTH separation
# ==================================

            sep_sol = ne_sep_sol.NE_SEP_SOL(nshot, self.tot)
            for lbl in ('LID_sol', 'LID_sep', 'IOC_sol', 'IOC_sep', 'NE_SOL', 'NE_SEP'):
                self.tth[lbl] = sep_sol.__dict__[lbl]
            for lbl in ('iflux', 'lid_flag', 'ioc_flag'):
                self.tth['loss_par'][lbl] = sep_sol.__dict__[lbl]

            self.tth['loss_par']['fac_sol']  = 0.1
            self.tth['loss_par']['fac_orb']  = 0
            self.tot['heat_par']['fac_ecrh'] = 1

# ============ 
# NBI losses
# ============ 

            self.tth['Te_av']   = np.zeros_like(ttot)
            self.tth['Ec']      = np.zeros_like(ttot)
            self.tth['tau_sp']  = np.zeros_like(ttot)
            self.tth['tau_sd']  = np.zeros_like(ttot)
            self.tth['P_i/P_e'] = np.zeros_like(ttot)
            tau_sd = np.zeros((nttot, nnbi))

# Auxiliary signals

            pmhd = self.tot['Wmhd']/np.maximum(self.tot['Vol'], 0.01)

            pnbi_max = np.max(self.tot['PNBI_TOT'])
            if pnbi_max > 0:
                print('\nCalculating slowing down time')

                indt = (pmhd > 100) & (ne_lav > 0.01)
                self.tth['Te_av'][indt]  = 0.2067*pmhd[indt]/ne_lav[indt]
                self.tth['Ec'][indt]     = 0.01865*self.tth['Te_av'][indt]
                self.tth['tau_sp'][indt] = np.maximum(7.376*1e-6*np.power(self.tth['Te_av'][indt], 1.5), 0)/ne_lav[indt]

                ind_src = (ni_eny > 0)
                for jsrc in range(nnbi):
                    tau_sd[indt, jsrc] = self.tth['tau_sp'][indt]/3*np.log(1 + np.power(ni_eny[jsrc]/self.tth['Ec'][indt], 1.5))
                tau_sd[:, ~ind_src] = 0.020
                for jt in range(nttot):
                    if (pmhd[jt] > 100) and (ne_lav[jt] > 0.01):
                        pie = np.zeros(nnbi, dtype=np.float32)
                        pie[ind_src] = self.tth['Ec'][jt]*3.0/(2.0*ni_eny[ind_src])
                        qie = 3.68*np.power(pie, 1.2)
                        fac = pniq[jt, :]/(1.0+qie)
                        prel_i = np.sum(qie[ind_src]*fac[ind_src])
                        prel_e = np.sum(fac[ind_src])
                        if prel_e > 5000:
                            self.tth['P_i/P_e'][jt] = prel_i/prel_e

                self.tth['tau_sd'][0] = 0.020
                for jt in range(1, nttot):
                    if self.tot['PNBI_TOT'][jt] > min_nbi:
                        p_src = pniq[jt, :]
                        ind_src = (p_src > min_nbi)
                        tmp = np.sum(p_src[ind_src]*tau_sd[jt, ind_src])
                        self.tth['tau_sd'][jt] = tmp/self.tot['PNBI_TOT'][jt]    
                    else:            # If NBI off, keep previous values of tau_sd
                        self.tth['tau_sd'][jt] = self.tth['tau_sd'][jt-1]
                    if self.tth['tau_sd'][jt] <= 0: # Keep inside loop, it influences next time step
                        self.tth['tau_sd'][jt] = 0.02

#Compute losses

            twfi = np.zeros_like(ttot)
            self.tth['SOL_LOSS'] = np.zeros_like(ttot)
            self.tth['RIP_LOSS'] = np.zeros_like(ttot)
            self.tth['CX_LOSS' ] = np.zeros_like(ttot)
            self.tth['ORB_LOSS'] = np.zeros_like(ttot)
            self.tot['SHINE_TH'] = np.zeros_like(ttot)
            
# Calculate NBI losses

            if pnbi_max > 0:
                print('\nCalculating NBI losses')
                if NBIpar == 'TRANSP':
                    print('Using TRANSP-based NBI parameterization')
                else:
                    print('Using FAFNER-based NBI parameterization')
                jsrc = 0
                Te_lav = self.tth['Te_av']*1.e-3 # keV
                for pn in ('INJ1', 'INJ2'):
                    for jnbi in range (0, n_nbi[pn]):
                        pnbi = pniq[:, jsrc]
                        if (dnis[pn]['UEXQ'][jnbi] <= 0) and ( np.max(pnbi) > min_nbi):
                            msg = 'Inconsistent NBI data'
                            msg_quit(nshot, msg)
                            return
                        if (np.max(pnbi) > min_nbi) and (dnis[pn]['UEXQ'][jnbi] > 0):
                            EdivA = dnis[pn]['UEXQ'][jnbi]/dnis[pn]['M']
                            self.tth['SOL_LOSS'] += pnbi*self.tth['NE_SOL']*1.e-19* \
                                     NBIlosses.solloss(nshot, pow_frac[pn], EdivA, jsrc)
                            self.tth['RIP_LOSS'] += NBIlosses.riploss(nshot, jsrc)
                            if NBIpar == 'FAFNER':
                                loss = NBIlosses.FAFPAR(nshot, pow_frac[pn], ne_lav, Te_lav, \
                                       self.tth['NE_SEP']*1e-19, tau_sd[:, jsrc], EdivA, jsrc)

                            if NBIpar == 'TRANSP':
                                loss = NBIlosses.TRAPAR(nshot, ipl_sign, Te_lav, 
                                       ne_lav, np.abs(self.tot['IpiFP'])*1.e-6, \
                                       EdivA, self.tot['peak'], self.tot['Zeff'], jsrc)

                            self.tth['CX_LOSS' ] += pnbi*loss.cx
                            self.tth['ORB_LOSS'] += pnbi*loss.orb
                            self.tot['SHINE_TH'] += pnbi*loss.st
                            twfi += pnbi*loss.wfi
                        jsrc += 1

            tot_loss = self.tot['SHINE_TH']
            for sig in ('SOL', 'ORB', 'CX', 'RIP'):
                tot_loss += self.tth['%s_LOSS' %sig]
 
# Modify Wfi retaining NBI switching

            self.tot['Wfi'] = wfi.wfi(twfi, self.tth['tau_sd'], dt)
            self.tot['Wth'] = self.tot['Wmhd'] - self.tot['Wfi']

            self.tot['P_TOT'] = self.tot['P_OH'] + self.tot['heat_par']['fac_icrh']*self.tot['PICR_TOT'] + \
                            self.tot['heat_par']['fac_ecrh']*self.tot['PECR_TOT'] + \
                            self.tot['PNBI_TOT'] - self.tot['SHINE_TH'] - self.tot['dWmhd/dt']
            tot_loss = np.minimum(tot_loss, self.tot['P_TOT'])
            tot_loss = np.maximum(tot_loss, 0.)
            self.tth['PNBI_NET'] = self.tot['PNBI_TOT'] - tot_loss
            self.tth['P_NET']    = self.tot['P_TOT']    - tot_loss
            self.tot['P_NET']    = self.tth['P_NET']    # Double; useful for scaling laws, deleted in the end
            self.tot['tau_tot']  = self.tot['Wmhd']/self.tot['P_TOT']
            self.tth['tau_th']   = self.tot['Wth' ]/self.tth['P_NET']
            self.tth['tau_hybr'] = self.tot['Wmhd']/self.tth['P_NET']

# Beta

            norm = 167.55*self.tot['ahor']/(self.tot['BTF']*self.tot['IpiFP']*self.tot['Vol'])
            self.tot['beta_N']   = norm*self.tot['Wmhd']
            self.tot['beta_Nth'] = norm*self.tot['Wth']

# Greenwald density

            self.tot['n/nGW'] = 1.e5*np.pi*self.tot['H-1_corr']*self.tot['ahor']**2/self.tot['IpiFP']
            ind = np.where(self.tot['ahor'] <= 0.05)
            self.tot['n/nGW'][ind] = 0

# Tau from scaling laws

            print('\nCalculating tau from scaling laws')
            tauscal = {}
            hfac = {}
            ion_mass = np.float32(ion_mass)
            for law in tth_laws[:-1]:
                print(law)
                parms = laws[law]
                tauscal[law] = np.zeros_like(ttot)
                tauscal[law] += parms['coeff']*np.power(ion_mass, parms['A'])

                for sig in ('IpiFP', 'Rgeo', 'ahor', 'kappa', 'karea', 'P_TOT', 'H-1_corr', 'BTF', 'P_NET'):
                    sig_fac = 1
                    if sig in ('IpiFP', 'P_TOT', 'P_NET'):
                        sig_fac = 1.e-6
                    ind_pos = (np.nan_to_num(self.tot[sig]) > 0)
                    tauscal[law][ind_pos] *= np.power(sig_fac*self.tot[sig][ind_pos], parms[sig])
                    tauscal[law][~ind_pos] = 0

                scal_mod = parms['scalmod']
                if scal_mod in self.tot.keys():
                    tau = self.tot[scal_mod]
                else:
                    tau = self.tth[scal_mod]

                ind_pos = (tauscal[law] > 0)
                hfac[law] = np.zeros_like(ttot)
                hfac[law][ind_pos] = tau[ind_pos]/tauscal[law][ind_pos]

# Kardaun scalings
# 1) Eq 2 in Kardaun O. et al., 
# Synopsis (IAEA 2006, Chengdu, China) | (IT/P1-10)
# http://efdasql.ipp.mpg.de/Igd/ITER/Conferences/IAEA2006/iaea06_syn35.pdf
# 2) Equation page 7 in P. T. Lang et al.
# IAEA 2012 | EX /P4-01
            self.tot['qcyl'] = 5*self.tot['Area']*self.tot['BTF']/(np.pi*self.tot['Rgeo']*1.e-6*self.tot['IpiFP'])
            self.tot['R/a']      = self.tot['Rgeo']/self.tot['ahor']
            self.tot['q95/qcyl'] = self.tot['q95']/self.tot['qcyl']

            kard = { \
                'IpiFP':1.4, 'A':0.2, 'BTF':0.12, 'P_NET':0.26, 'Rgeo':1.17, \
                'R/a':(0.51, -0.90), 'karea':0.37, 'q95/qcyl':0.77, 'n/nGW':(0.11, -0.22) }
            kard_lang = { \
                'IpiFP':1.4, 'A':0.2, 'BTF':0.12, 'P_NET':0.26,    'Rgeo':1.17, \
                'R/a':(2.546, -0.9), 'karea':0.37, 'q95/qcyl':0.77, 'n/nGW':(0.0322, -0.22) }

# q95 = 3
# Vol = 831
# Area = 21.33
# karea = Area/(pi * a**2) = 1.698
# qcyl = 5 Bt*area/(Ip pi Rgeo) = 5 Bt * vol/(2 Ip pi**2 Rgeo**2) = 0.2533 *Bt*Vol/(Ip Rgeo**2) =   1.935

            iter_ref = {'IpiFP':15e6, 'A':2.5, 'BTF':5.3, 'P_NET':87.7e6, 'Rgeo':6.2, 'R/a':3.1, \
                        'karea':1.698, 'q95/qcyl':1.550, 'n/nGW':0.838 }

            lnWK  = np.zeros_like(ttot)
            logpn = np.zeros_like(ttot)
            lnWK += 5.6348 #ln 280
            lnWK += kard['A']*np.log(ion_mass/iter_ref['A'])
            for pn in('IpiFP', 'BTF', 'Rgeo', 'karea', 'q95/qcyl', 'P_NET'):
                loc_dic = np.nan_to_num(self.tot[pn])
                ind_pos = (loc_dic > 0)
                logpn[ind_pos] = np.log(loc_dic[ind_pos]/iter_ref[pn])
                logpn[~ind_pos] = 1.e-5
                lnWK += kard[pn]*logpn

            for pn in('n/nGW', 'R/a'):
                ind_pos = (self.tot[pn] > 0)
                logpn[ind_pos] = np.log(self.tot[pn][ind_pos]/iter_ref[pn])
                logpn[~ind_pos] = 1.e-5
                lnWK += kard[pn][0]*logpn + kard[pn][1]*logpn**2

            law = 'KARDAUN, IAEA 2006, IT/P10'
            tauscal[law] = 1.e6*np.exp(lnWK)/self.tth['P_NET']
            tmp2 = self.tth['tau_th']/tauscal[law]
            tmp2[np.isinf(tmp2)] = 0
            hfac[law] = np.nan_to_num(tmp2)

            self.tth['Scalings'] = np.zeros((nttot, len(tth_laws)), dtype=np.float32)
            self.tth['H/L-facs'] = np.zeros((nttot, len(tth_laws)), dtype=np.float32)

            for pn in parlaws.keys():
                self.tth['scal_par'][pn] = np.zeros(len(tth_laws), dtype=np.float32)

# H factors

            for jlaw, law in enumerate(tth_laws):
                self.tth['H/L-facs'][:, jlaw] = hfac[law]
                self.tth['Scalings'][:, jlaw] = tauscal[law]
                for pn, base in parlaws.items():
                    self.tth['scal_par'][pn][jlaw] = laws[law][base]

            self.tth['scal_par']['descript'] = tth_laws

# L->H scaling
# Y. Martin, Journal of Physics: Conference Series 123 (2008) 012033, eq 2-3
            lh_laws = { \
                'ITPA Threshold Power Scaling 2007 (Plasma Surface S)': \
                {'prefac':0.0488, 'Rgeo':0, 'ahor':0, 'S':0.941, '<ne>_H-1':0.717, 'Btf':0.803}, \
                'ITPA Threshold Power Scaling 2007 (Rgeo and ahor)': \
                {'prefac':2.150, 'Rgeo':0.999, 'ahor':0.975, 'S':0, '<ne>_H-1':0.782, 'Btf':0.772} \
            }
            lhlaws = ('ITPA Threshold Power Scaling 2007 (Plasma Surface S)', \
                      'ITPA Threshold Power Scaling 2007 (Rgeo and ahor)')

            nlh = len(lhlaws)
            self.tth['L2H_SCAL'] = np.zeros((nttot, nlh), dtype=np.float32)
            self.tth['L2H_facs'] = np.zeros((nttot, nlh), dtype=np.float32)
            self.tth['L2H_par'] = {}
            l2hscal = {}
            l2hfac  = {}

            for pn in ('prefac', 'Rgeo', 'ahor', 'S', '<ne>_H-1', 'Btf'):
                self.tth['L2H_par'][pn] = np.zeros(nlh, dtype=np.float32)

            for law, par in lh_laws.items():
                l2hscal[law] = np.zeros_like(ttot)
                l2hfac[law]  = np.zeros_like(ttot)
                tmp = np.ones(nttot)
                for sig in ('Rgeo', 'ahor', '<ne>_H-1', 'Btf'):
                    if sig == '<ne>_H-1':
                        sig_fac = 1.e-1
                    else:
                        sig_fac = 1
                    arr = self.tot[parlaws[sig]]
                    ind_pos = (par > 0)
                    tmp[ind_pos] *= np.power(sig_fac*arr[ind_pos], par[sig])
                    tmp[~ind_pos] = 0
                ind_pos = (self.tot['Surface'] > 0)
                tmp[ind_pos] *= 1e6*par['prefac']*np.power(self.tot['Surface'][ind_pos], par['S'])
                tmp[~ind_pos] = 0
                l2hscal[law] = tmp
                ind = np.where(tmp > 0)
                l2hfac[law][ind] = self.tth['P_NET'][ind]/tmp[ind]

            for jlaw, law in enumerate(lhlaws):
                self.tth['L2H_facs'][:, jlaw] = l2hfac[law]
                self.tth['L2H_SCAL'][:, jlaw] = l2hscal[law]
                for pn, val in lh_laws[law].items():
                    self.tth['L2H_par'][pn][jlaw] = val
            self.tth['L2H_par']['descript'] = lhlaws

# Q/(Q+5) figure of merit

            Ggw_laws = { \
                'ITERH-98P(y, th, 2), ELMy, Wth/Pnet': \
                {'prefac':35.5, 'G_H_exp':3.23, 'G_bn_exp':-1.23, 'G_q_exp':-3.1}, \
                'DS03 or ESGB, McDonald, ppcf 46 (2004) A215 equation 11': \
                {'prefac':22.8, 'G_H_exp':2.22, 'G_bn_exp':-0.22, 'G_q_exp':-2.71}, \
                'CORDEY05, NF 45 (2005) 1078, equation 9': \
                {'prefac':6.07, 'G_H_exp':1.82, 'G_bn_exp':0.18, 'G_q_exp':-2.22} \
            }

            nGgw = len(Ggw_laws)
            self.tth['G_gw'] = np.zeros((nttot, nGgw), dtype=np.float32)
            self.tth['G_gw_par'] = {}
            for pn in ('prefac', 'G_H_exp', 'G_bn_exp', 'G_q_exp'):
                self.tth['G_gw_par'][pn] = np.zeros(nGgw, dtype=np.float32)

            jlaw = 0
            for law, par in Ggw_laws.items():
                ind_pos = (hfac[law] > 0) & (self.tot['beta_Nth'] > 0) & (self.tot['q95'] > 0)
                self.tth['G_gw'][ind_pos, jlaw] = par['prefac'] * \
                    np.power(hfac[law][ind_pos], par['G_H_exp']) * \
                    np.power(self.tot['beta_Nth'][ind_pos], par['G_bn_exp']) * \
                    np.power(self.tot['q95'][ind_pos], par['G_q_exp'])
                for pn, val in par.items():
                    self.tth['G_gw_par'][pn][jlaw] = val
                jlaw += 1
            self.tth['G_gw_par']['descript'] = (tth_laws[7], tth_laws[8], tth_laws[9])

# Final settings

            self.tth['TOT_file']['expr'] = exp_write
            self.tth['TOT_file']['edition'] = -1
            if sf.Open('TOT', nshot):
                print('TOT edition used for TTH: %d' %sf.edition)
                self.tth['TOT_file']['edition'] = sf.edition

            self.tot['H-1_corr'] *= 1e19
            self.tot['H-4_corr'] *= 1e19

# Delete some entries

            for sig in ('Prad', 'P_NET'):
                if sig in self.tot.keys():
                    del self.tot[sig]

        except ZeroDivisionError:
            print('\nDivision by zero\n')
            raise
        except NameError:
            print('\nVariable not defined\n')
            raise
        except TypeError as TE:
            print('\nType error\n')
#            return
            raise
        except ValueError:
            print('\nValue error\n')
            raise
        except IOError:
            print('\nCannot open\n')
            raise
        except:
            print('\nUnexpected error\n')
            raise
