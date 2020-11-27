import neut_flux
import numpy as np
import dd_20180130

sf  = dd_20180130.shotfile()


class NE_SEP_SOL:


    def __init__(self, nshot, tot_d):
# ================================== 
# Estimate ne at SOL, sep
# only reason for TOT-TTH separation
# ==================================

# LID
        diag = 'LID'
        dlid = {}
        lid_flag = 1
        if sf.Open(diag, nshot):
            lid_flag = 0
            tlid = sf.GetTimebase('ne_sol')
            if type(tlid) == type(np.float64(0)):
                lid_flag = 1
            for sig in ('ne_sol', 'ne_sep'):
                dlid[sig] = sf.GetSignal(sig, cal=True)
            if dlid['ne_sol'] is None:
                print('LID:%s not found, discarding LID data' %sig)
                lid_flag = 1
            elif np.isnan(dlid['ne_sol']).any():
                print('LID:%s has a NaN, discarding LID data' %sig)
                lid_flag = 1
            elif np.max(dlid['ne_sol']) <= 0:
                print('LID:%s has no positive values, discarding LID data' %sig)
                lid_flag = 1
            sf.Close()

# DLP
        diag = 'DLP'
        ddlp = {}
        if sf.Open(diag, nshot):
            tdlp = sf.GetTimebase('ne_sol')
            if tdlp is not None:
                for sig in ('ne_sol', 'ne_sep'):
                    ddlp[sig] = sf.GetSignal(sig, cal=True)
                sf.Close()
                if (lid_flag == 1): #LID not available
                    tlid = tdlp
                    dlid = ddlp
                else:
                    for sig in ('ne_sol', 'ne_sep'):
                        (ind_lid_zero, ) = np.where(dlid[sig] == 0)
# Where ne_sol or ne_sep is zero, replace by DLP value at the same time
                        for jt in ind_lid_zero:
                            print('tlid jt', tlid[jt])
                            if (tlid[jt] == tdlp.any()):
                                jtim = np.min(np.where(tdlp == tlid[jt]))
                                dlid[sig][jt] = ddlp[sig][jtim]

# IOC

        neut = neut_flux.IOC_FLUX(nshot)

        self.iflux    = neut.iflux
        self.lid_flag = lid_flag
        self.ioc_flag = neut.ioc_flag

        ntlid = 0
        if lid_flag == 0:
            ntlid = len(tlid)
        dtlid_max = 0.15
        if (neut.ioc_flag == 1) and (lid_flag == 1):
            tlid = np.array([0, 10])

# If LID and IOC, interp LID onto tioc
# set netflgc[0] = 0 where tioc is outside tlid or where delta_tlid > 0.15 s
# Old divertor

        ttot = tot_d['time']

        self.LID_sol = np.zeros_like(ttot)
        self.LID_sep = np.zeros_like(ttot)
        self.IOC_sol = np.zeros_like(ttot)
        self.IOC_sep = np.zeros_like(ttot)
        self.NE_SOL  = np.zeros_like(ttot)
        self.NE_SEP  = np.zeros_like(ttot)

        if neut.ioc_flag == 0: # IOC ok, first choice
            diag = 'BPD'
            sig = 'Prad'
            if sf.Open(diag, nshot):
                tbpd = sf.GetTimebase(sig)
                prad = sf.GetSignal(sig)
                if prad is None:
                    return
                else:
                    prad = np.nan_to_num(prad)
                sf.Close()
            else:
                return

            dhelp = {}

            dhelp['Prad'] = np.interp(neut.tioc, tbpd, prad)
            for sig in ('IpiFP', 'H-1_corr', 'q95', 'PNBI_TOT', 'PECR_TOT', 'PICR_TOT', 'u_loop'):
                dhelp[sig] = np.interp(neut.tioc, ttot, tot_d[sig])
            dhelp['H-1_corr'] *= 1e19
            if nshot < 8646:
# Fit after Joe Schweinzer (EPS 96), Divertor I
                phelp = dhelp['PNBI_TOT'] + dhelp['PICR_TOT'] + dhelp['PECR_TOT']
                phelp *= 1e-6
                nesol_ioc = np.zeros_like(phelp)
                nesep_ioc = np.zeros_like(phelp)
                ind_sol = (dhelp['Prad'] > 0) & (neut.neut_flux > 0)
                ind_sep = ind_sol & (phelp > 0)
                nesol_ioc[ind_sol] = 7e5*np.power(neut.neut_flux[ind_sol], 0.51) * \
                              np.power(np.abs(dhelp['q95'][ind_sol]), 0.4) * \
                              np.power(dhelp['Prad'][ind_sol], 0.24)
                nesep_ioc[ind_sep] = np.exp(21)*np.power(neut.neut_flux[ind_sep], 0.32) * \
                                     np.power(np.abs(dhelp['IpiFP'][ind_sep]), 0.47) * \
                                     np.power(dhelp['Prad'][ind_sep]*1e-6, 0.172) * \
                                     np.power(phelp, 0.12)
            else:
# Fit after Joe Schweinzer, February 1999 for Divertor II
                phelp = 0.9*dhelp['PNBI_TOT'] + 0.8*dhelp['PICR_TOT'] + dhelp['PECR_TOT'] + \
                            dhelp['IpiFP']*dhelp['u_loop'] - dhelp['Prad']
                nesol_ioc = np.zeros_like(phelp)
                nesep_ioc = np.zeros_like(phelp)
                ind_pos = (phelp > 0) & (dhelp['H-1_corr'] > 0) & (neut.neut_flux > 0)

                nesol_ioc[ind_pos] = np.exp(-0.78) * np.power(neut.neut_flux[ind_pos], 0.42) * \
                                          np.power(np.abs(dhelp['q95'][ind_pos]), 0.57) * \
                                          np.power(phelp[ind_pos], -0.05) * \
                                          np.power(dhelp['H-1_corr'][ind_pos], 0.49)
                nesep_ioc[ind_pos] = np.exp(1.32) * np.power(neut.neut_flux[ind_pos], 0.16) * \
                                         np.power(phelp[ind_pos], 0.03) * \
                                         np.power(dhelp['H-1_corr'][ind_pos], 0.76)

# Set inf, NaN to zero
                self.IOC_sol = np.interp(ttot, neut.tioc, nesol_ioc)
                self.IOC_sep = np.interp(ttot, neut.tioc, nesep_ioc)

        if lid_flag == 0:
            indpos = (dlid['ne_sol'] > 0)
            nesol = dlid['ne_sol'][indpos]
            tlid_red = tlid[indpos]
            self.LID_sol = np.interp(ttot, tlid_red, nesol)
            indpos = (dlid['ne_sep'] > 0)
            nesep = dlid['ne_sep'][indpos]
            tlid_red = tlid[indpos]
            self.LID_sep = np.interp(ttot, tlid_red, nesep)

        if neut.ioc_flag == 0:
            self.NE_SOL = self.IOC_sol
            self.NE_SEP = self.IOC_sep
        elif lid_flag == 0:
            print('Using LID for NE_SOL')
            self.NE_SOL = self.LID_sol
            self.NE_SEP = self.LID_sep
        else:
            print('Ill defined LID')
