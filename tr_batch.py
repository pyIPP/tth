import os, sys
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')
sys.path.append('/afs/ipp/home/g/git/trview')
import numpy as np
from multiprocessing import Pool, cpu_count
import read_prof, time_slice, main_spec, fit_sep
import get_sig, get_ech, get_nbi, get_icrf
from vprof import fit1d
from sf2equ_20200525 import EQU

pre_list = ['E', 'I', 'N', 'V']

fit_d = {'I': 'rec_spl', 'E': 'rec_spl', 'N': 'rec_spl', 'V': 'rec_spl'}

prio = { \
    'I': ['CEZ:Ti_c', 'CEZ:Ti', 'CHZ:Ti_c', 'CHZ:Ti'], \
    'V': ['CEZ:vrot', 'CHZ:vrot'], \
    #'E': ['IDA:Te', 'VTA:Te_c', 'CEC:Trad-A', 'VTA:Te', 'VTN:Te'], \
    'E': ['IDA:Te', 'VTA:Te_c', 'VTA:Te', 'VTN:Te', 'CEC:Trad-A'], \
    'N': ['IDA:ne', 'VTA:Ne_c', 'VTA:Ne', 'VTN:Ne', 'DLP:ne', 'DPR:ne', 'YAG:ne'] \
}

tr_dir = '%s/tr_client/AUGD' %os.getenv('HOME')


class PROF:

    pass


class GPAR:


    def __init__(self, nshot, tbeg, tend):

        self.shot    = nshot
        self.tbeg    = tbeg
        self.tend    = tend
        self.eqdia   = ['EQI', 'EQH', 'FPP']
        self.nmom    = 6
        self.eq_exp  = 'AUGD'
        self.eq_ed   = 0
        self.nrho    = 51
        self.dt_prof = 0.01
        self.angf    = True
        self.no_elm  = True
        self.fit_tol = 0.3
        self.sigmas  = 1.


class INPUT:


    def __init__(self, gpar, code='rb'):

        self.gpar = gpar
        self.get_data(code=code)


    def get_data(self, code='rb'):

        gpar = self.gpar
        self.nbi = get_nbi.NBI  (gpar.shot, tbeg=gpar.tbeg, tend=gpar.tend)
        if code in ('tr', 'as'):
            self.ech = get_ech.gyros(gpar.shot, tbeg=gpar.tbeg, tend=gpar.tend)
            self.icr = get_icrf.ICRF(gpar.shot, tbeg=gpar.tbeg, tend=gpar.tend)

        for diag in gpar.eqdia:

            try:

                self.equ = EQU(gpar.shot, tbeg=gpar.tbeg, tend=gpar.tend, diag=diag, exp=gpar.eq_exp, ed=gpar.eq_ed)
                if code in ('tr', 'as'):
                    self.equ.geofit = fit_sep.geofit(self.equ, nmom=gpar.nmom, n_the=gpar.n_the)
                break
                
            except:
                print('Equil. diag %s not found' %diag)
                continue
                
        self.sig = get_sig.READ(self.equ)
        self.load_2d_profiles(gpar)


    def load_2d_profiles(self, gpar):

        diags, sigs, pres, err_sigs = \
            np.genfromtxt('/afs/ipp/home/g/git/trview/diags.dat', dtype='str', skip_header=2, unpack=True, usecols=(0, 1, 2, 6))

        raw_data = {}
        x_fit = np.linspace(0, 1, gpar.nrho)

# Fetch exp data

        raw_data = {}
        sfl_d = {}
        for prefix in pre_list:
            sfl_d[prefix] = []
            for diag_sigs in prio[prefix]: #['CEZ:Ti_c CMZ:Ti_c', 'CHZ:Ti_c']
                for diag_sig in diag_sigs.split(): #'CEZ:Ti_c CMZ:Ti_c'
                    diag_sig = diag_sig.strip()
                    sflist = 'AUGD:%s:0' %diag_sig
                    rdd = read_prof.dd2read(sflist, self.equ, gpar)
                    if rdd is not None:
# Replace edition 0 with actual editoin number
                        sf_list = sflist.split(':')
                        sfl = ':'.join(sf_list[:-1])
                        sfl += ':%d' %rdd['ed']
                        raw_data[sfl] = rdd
                        sfl_d[prefix].append(sfl)
                if len(sfl_d[prefix]) > 0: # priority, go to next only if first fails
                    break

        if len(sfl_d['E']) == 0 or len(sfl_d['N']) == 0:
            raise Exception('Missing Te or ne profile, no run possible')

# Create wish time grid

        tbeg = 10.
        tend = 0.
        for rdd in raw_data.values():
            tbeg = min(tbeg, rdd['tgrid'][0])
            tend = max(tend, rdd['tgrid'][-1])
        tbeg = max(tbeg, gpar.tbeg)
        tend = min(tend, gpar.tend)

        nt = int((tend - tbeg)/ gpar.dt_prof)
        tgrid = tbeg + gpar.dt_prof*np.arange(nt)

# Profile object: start with parameter settings

        self.profiles = gpar
        n_spl = self.profiles.nrho
        self.profiles.time = tgrid
        self.profiles.rho  = np.linspace(0, 1, n_spl)
        self.profiles.prof = {}

# slice all kinetic profiles

        for prefix in pre_list:
            red_data = {}
            for sflist in sfl_d[prefix]:
                red_data[sflist] = {}
                x_red, y_red, err_red = time_slice.map2tgrid( \
                raw_data[sflist], tgrid, nshot=gpar.shot, noELMs=gpar.no_elm)
                red_data[sflist]['rhop']     = x_red
                red_data[sflist]['data']     = y_red
                red_data[sflist]['data_err'] = err_red
                red_data[sflist]['tgrid']    = tgrid
                red_data[sflist]['unit']     = raw_data[sflist]['unit']

# Combine diagnostics

            xexp = None
            for sflist in sfl_d[prefix]:
                if xexp is None:
                    xexp = np.array(red_data[sflist]['rhop'])
                    yexp = np.array(red_data[sflist]['data'])
                    yerr = np.array(red_data[sflist]['data_err'])
                else:
                    xexp = np.append(xexp, red_data[sflist]['rhop']    , axis=1)
                    yexp = np.append(yexp, red_data[sflist]['data']    , axis=1)
                    yerr = np.append(yerr, red_data[sflist]['data_err'], axis=1)

            print('prefix:', prefix)

            if xexp is not None:

                print(xexp.shape)
                ind_sort = np.argsort(xexp, axis=1)
                for jt in range(xexp.shape[0]):
                    xexp[jt, :] = xexp[jt, ind_sort[jt]]
                    yexp[jt, :] = yexp[jt, ind_sort[jt]]
                    yerr[jt, :] = yerr[jt, ind_sort[jt]]
# Fit


                if fit_d[prefix] == 'rec_spl':
                    timeout_pool = 10
                else:
                    timeout_pool = 200

                print('Fitting', sfl_d[prefix])
                if prefix in ('E', 'N'):
                    fit_tol = 0.
                else:
                    fit_tol = gpar.fit_tol

                try:
                    pool = Pool(cpu_count())
                    out = pool.map_async(fit1d, [(xexp[jt], yexp[jt], yerr[jt], \
                      x_fit, fit_tol, gpar.sigmas, fit_d[prefix]) \
                      for jt in range(nt)]).get(timeout_pool)
                    pool.close()
                    pool.join()
        
                    probj = PROF()
                    probj.fit_rhop = np.array(out)[:, :n_spl]
                    probj.sflists  = sorted(red_data.keys())
                    probj.fit_method = fit_d[prefix]
                    probj.tnot = []
                    probj.units = red_data[probj.sflists[0]]['unit']
                    self.profiles.prof[prefix] = probj
                except:
                    print('Problems fitting prefix %s' %prefix)

        if 'I' not in self.profiles.prof:

            print('\n No experimental data for Ti, setting Ti = Te \n')

            self.profiles.prof['I'] = self.profiles.prof['E']

        if 'V' not in self.profiles.prof:

            print('\n\n Setting Vrot to zero \n\n')
            probj = PROF()
            probj.fit_rhop = 0*np.array(out)[:, :n_spl]
            probj.sflists  = sorted(red_data.keys())
            probj.fit_method = fit_d['V']
            probj.tnot = []
            probj.units = 'rad/s'
            self.profiles.prof['V'] = probj

        self.ext_d = {}
        for prefix, sfls in sfl_d.items():
            self.ext_d[prefix] = ''
            for sfl in sfls:
                diag = sfl.split(':')[1]
                self.ext_d[prefix] += diag

        ion_spec = main_spec.spec_jou_sf(gpar.shot)
        print('Main ion species is %s' %ion_spec)
        self.ion_d = {'Aimp': 12, 'Zimp': 6}
        if ion_spec == 'D':
            self.ion_d['Amain'] = 2
            self.ion_d['Zmain'] = 1
        elif ion_spec == 'H':
            self.ion_d['Amain'] = 1 
            self.ion_d['Zmain'] = 1
        elif ion_spec == 'He':
            self.ion_d['Amain'] = 4 
            self.ion_d['Zmain'] = 2 
