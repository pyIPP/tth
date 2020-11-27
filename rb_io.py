import os, sys, shutil, datetime
import read_beam
import numpy as np
from scipy.io import netcdf
sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib')
from sfarr import SFARR


if os.getenv('TRVIEW_RBHOME') is None:
    rbhome = '/toks/work/%s/rabbit' %os.getenv('USER')
else:
    rbhome = os.getenv('TRVIEW_RBHOME')


def rb2sf(rb_in, rb_out):

    import ww_20180130 as ww
    import sfh_20200703 as sfh

# Paths

    homdir = os.getenv('HOME')
    os.system('mkdir -p %s/shotfiles/RAB' %homdir)
    source = '/afs/ipp/home/g/git/python/rbview/RAB00000.sfh.temp'
    sfhdir = '%s/python/rbview' %homdir
    os.system('mkdir -p %s' %sfhdir)

    fsfh = '%s/RAB00000.sfh' %sfhdir

    try:
        shutil.copy2(source, fsfh)
    except:
        print('Unable to copy file %s to %s' %(source, fsfh))
        raise Exception

    sfh_dic = sfh.sfhmod(rb_in, rb_out, nml='', fsfh=fsfh)
    sfh_dic['torqjxb'] = np.sum(sfh_dic['torqjxb'], axis=2)
    ww.write_sf(rb_out.shot, sfh_dic, sfhdir, 'RAB', exp=os.getenv('USER'))


def rb2cdf(rb_in, rb_out):

    f_cdf = '%s/%drb.cdf' %(rbhome, rb_out.shot)
    f = netcdf.netcdf_file(f_cdf, 'w', mmap=False)

    f.history = "Created " + datetime.datetime.today().strftime("%d/%m/%y")

    nt, nrho_out, nnbi, nv = rb_out.bdepQ.shape
    nrho_in = len(rb_in.rho_in)
    nrho_eq = rb_in.psi1d.shape[1]

    f.createDimension('time'   , nt)
    f.createDimension('rho_out', nrho_out)
    f.createDimension('rho_in' , nrho_in)
    f.createDimension('rho_eq' , nrho_eq)
    f.createDimension('spec'   , nv)
    f.createDimension('n_nbi'  , nnbi)

    for sig, sf_arr in rb_out.__dict__.items():
        if sig in ('shot', ):
            continue
#        print(sig, sf_arr.shape, sf_arr.dims)
        tmp = f.createVariable(sig, np.float32, sf_arr.dims)
        tmp[:]    = sf_arr
        tmp.units = sf_arr.units
        tmp.long_name = sf_arr.long_name

    for sig in ('rho_in', 'te', 'ti', 'ne', 'wtor', 'zeff', 'pnbi', \
                'psi_sep', 'psi_axis', \
                'rho1d', 'psi1d', 'q', 'F', 'vol', 'area'):
        sf_arr = rb_in.__dict__[sig]
        tmp = f.createVariable(sig, np.float32, sf_arr.dims)
        tmp[:]    = sf_arr
        tmp.units = sf_arr.units
        tmp.long_name = sf_arr.long_name

    f.close()
    print('Stored %s' %f_cdf)
    return f_cdf


class RABBIT_IN:


    def __init__(self, nshot):


        nshot = int(nshot)
        self.shot = nshot
        self.rbpath = '%s/%d' %(rbhome, nshot)

        self.beam_data()
        if not hasattr(self, 'start_pos'):
            print('Input not found for shot %d' %nshot)
            return
        self.options()
        flim = self.opt['physics']['limiter_file']
        self.Rlim, self.Zlim = np.loadtxt(flim, skiprows=1, unpack=True)
        self.Rlim = SFARR(self.Rlim)
        self.Zlim = SFARR(self.Zlim)
        self.Rlim.units = '(m)'
        self.Zlim.units = '(m)'
        self.Rlim.dims = ('vessel-polygon', )
        self.Zlim.dims = ('vessel-polygon', )

        self.signals()
        self.equ(len(self.time))


    def options(self):

        fopt = '%s/options.nml' %self.rbpath
        self.opt = {}
        if not os.path.isfile(fopt):
            print('File %s not found' %fopt)
            return
        with open(fopt, 'r') as f:
            lines = f.readlines()
        for jlin in range(len(lines)):
            sline = lines[jlin].strip()
            if sline == '':
                continue
            if sline[0] == '&':
                key = sline[1:]
                self.opt[key] = {}
                while True:
                    jlin += 1
                    sline = lines[jlin].strip()
                    if sline[0] == '/':
                        break
                    key1, val1 = sline.split('=')
                    key1 = key1.strip()
                    try:
                        val = int(val1)
                    except:
                        try: 
                            val = float(val1)
                        except:
                            val = val1.strip().replace("'", '')
                    self.opt[key][key1] = val


    def beam_data(self):

        fbeam = '%s/beams.dat' %self.rbpath
        if not os.path.isfile(fbeam):
            print('File %s not found' %fbeam)
            return
        with open(fbeam, 'r') as f:
            lines = f.readlines()
        for jlin in range(len(lines)):
            line = lines[jlin].strip()
            if 'sources' in line:
                jlin += 1
                nnbi = int(lines[jlin])
            if ' nv' in line:
                jlin += 1
                nv = int(lines[jlin])
#
            if 'start pos' in line:
                self.start_pos = SFARR(np.zeros((nnbi, 3), dtype=np.float32))
                for jnb in range(nnbi):
                    self.start_pos[jnb] = lines[jlin + jnb + 1].split()
                jlin += nnbi

            if 'unit vector' in line:
                self.unit_vec = SFARR(np.zeros((nnbi, 3), dtype=np.float32))
                for jnb in range(nnbi):
                    self.unit_vec[jnb] = np.array([lines[jlin + jnb + 1].split()])
                jlin += nnbi

            if 'polynomial' in line:
                self.beam_width = SFARR(np.zeros((nnbi, 3), dtype=np.float32))
                for jnb in range(nnbi):
                    self.beam_width[jnb] = np.array([lines[jlin + jnb + 1].split()])
                jlin += nnbi
#
            if 'Injection energy' in line:
                self.einj = SFARR(np.zeros(nnbi, dtype=np.float32))
                for jnb in range(nnbi):
                    self.einj[jnb] = float(lines[jlin + jnb + 1])
                jlin += nnbi
            if 'A beam' in line:
                self.Abeam = SFARR(np.zeros(nnbi, dtype=np.float32))
                for jnb in range(nnbi):
                    self.Abeam[jnb] = float(lines[jlin + jnb + 1])
                jlin += nnbi
#
            if 'fraction' in line:
                self.part_frac = SFARR(np.zeros((nnbi, nv), dtype=np.float32))
                for jnb in range(nnbi):
                    self.part_frac[jnb] = np.array([lines[jlin + jnb + 1].split()])
                jlin += nnbi

        self.start_pos.units  = '(m)'
        self.unit_vec.units   = ''
        self.beam_width.units = ''
        self.einj.units       = '(eV)'
        self.Abeam.units      = '(AMU)'
        self.part_frac.units  = ''

        self.start_pos.long_name  = 'Position of NBI source in real space'
        self.unit_vec.long_name   = 'Versor of direction of central beam'
        self.beam_width.long_name = 'Polynomial coefficients for beam-width'
        self.einj.long_name       = 'Injection energy'
        self.Abeam.long_name      = 'Mass number of NBI ions'
        self.part_frac.long_name  = 'Particle fraction of full/half/third energy'

        self.start_pos.dims  = '(n_nbi, space)'
        self.unit_vec.dims   = '(n_nbi, space)'
        self.beam_width.dims = '(n_nbi, coeff)'
        self.einj.dims       = '(n_nbi, )'
        self.Abeam.dims      = '(n_nbi, )'
        self.part_frac.dims  = '(n_nbi, spec)'


    def signals(self):

        fsig = '%s/timetraces.dat' %self.rbpath
        with open(fsig, 'r') as f:
            lines = f.readlines()
        nlin = len(lines)
        nt = int(lines[0])
        nrho = int(lines[1])
        rho_lbl = lines[2].strip()
        data_block = []
        for jlin in range(3, nlin):
            data_block += lines[jlin].split()
        data_block = np.array(data_block, dtype=np.float32)
        jbeg = 0
        jend = nt
        self.time = SFARR(data_block[jbeg: jend])
        jbeg = nt
        jend = nt + nrho
        self.rho_in = SFARR(data_block[jbeg: jend])
        jbeg = jend
        jend = jbeg + nt*nrho
        self.te = SFARR(data_block[jbeg: jend]).reshape((nrho, nt)).T
        jbeg = jend
        jend = jbeg + nt*nrho
        self.ti = SFARR(data_block[jbeg: jend]).reshape((nrho, nt)).T
        jbeg = jend
        jend = jbeg + nt*nrho
        self.ne = SFARR(data_block[jbeg: jend]).reshape((nrho, nt)).T
        jbeg = jend
        jend = jbeg + nt*nrho
        self.wtor = SFARR(data_block[jbeg: jend]).reshape((nrho, nt)).T
        jbeg = jend
        jend = jbeg + nt*nrho
        self.zeff = SFARR(data_block[jbeg: jend]).reshape((nrho, nt)).T
        jbeg = jend
        self.pnbi = SFARR(data_block[jbeg: ])
        nnbi = len(self.pnbi)//nt
        self.pnbi = self.pnbi.reshape((nnbi, nt)).T

        self.time.units = '(s)'
        self.rho_in.units  = ''
        self.te.units   = '(keV)'
        self.ti.units   = '(keV)'
        self.ne.units   = '(1/cm**3)'
        self.wtor.units   = '(rad/s)'
        self.zeff.units  = ''
        self.pnbi.units = '(W)'

        self.time.long_name   = 'Time grid'
        self.rho_in.long_name = rho_lbl
        self.te.long_name   = 'Electron temperature'
        self.ti.long_name   = 'Ion temperature'
        self.ne.long_name   = 'Electron density'
        self.wtor.long_name   = 'Toroidal angular frequency'
        self.zeff.long_name  = 'Effective charge'
        self.pnbi.long_name = r'$P_{NBI}$'
        self.time.dims   = ('time', )
        self.rho_in.dims = ('rho_in', )
        self.pnbi.dims   = ('time', 'n_nbi')
        for attr in ('te', 'ne', 'ti', 'wtor', 'zeff'):
            self.__dict__[attr].dims   = ('time', 'rho_in')


    def equ(self, nt):

        self.psi_sep  = SFARR(np.zeros(nt))
        self.psi_axis = SFARR(np.zeros(nt))
        self.sgn_ip   = SFARR(np.zeros(nt))
        self.Rmag     = SFARR(np.zeros(nt))
        self.Zmag     = SFARR(np.zeros(nt))

        for jt in range(nt):

            fequ = '%s/equ/equ_%d.dat' %(self.rbpath, jt+1)
            if not os.path.isfile(fequ):
                break
            with open(fequ, 'r') as f:
                lines = f.readlines()

            nlin = len(lines)
            nr = int(lines[0])
            nz = int(lines[1])

            data_block = []
            for jlin in range(2, nlin):
                data_block += lines[jlin].split()
            data_block = np.array(data_block, dtype=np.float32)

            if jt == 0:
                self.Rmesh = SFARR(np.zeros(nr))
                self.Zmesh = SFARR(np.zeros(nz))
                self.psi2d = SFARR(np.zeros((nt, nz, nr)))
                self.rho2d = SFARR(np.zeros((nt, nz, nr)))

                self.Rmesh = SFARR(data_block[0: nr])
                self.Zmesh = SFARR(data_block[nr: nr + nz])

            jbeg = nr + nz
            jend = jbeg + nr*nz
            self.psi2d[jt] = data_block[jbeg: jend].reshape(nz, nr)
            jbeg = jend
            jend = jbeg + nr*nz
            self.rho2d[jt] = data_block[jbeg: jend].reshape(nz, nr)

            nrho = int(data_block[jend])
            if jt == 0:
                self.psi1d = SFARR(np.zeros((nt, nrho)))
                self.vol   = SFARR(np.zeros((nt, nrho)))
                self.area  = SFARR(np.zeros((nt, nrho)))
                self.rho1d = SFARR(np.zeros((nt, nrho)))
                self.q     = SFARR(np.zeros((nt, nrho)))
                self.F     = SFARR(np.zeros((nt, nrho)))
            jbeg = jend + 1
            jend = jbeg + nrho
            self.psi1d[jt] = data_block[jbeg: jend]
            jbeg = jend
            jend = jbeg + nrho
            self.vol[jt] = data_block[jbeg: jend]
            jbeg = jend
            jend = jbeg + nrho
            self.area[jt] = data_block[jbeg: jend]
            jbeg = jend
            jend = jbeg + nrho
            self.rho1d[jt] = data_block[jbeg: jend]
            jbeg = jend
            jend = jbeg + nrho
            self.q[jt] = data_block[jbeg: jend]
            jbeg = jend
            jend = jbeg + nrho
            self.F[jt] = data_block[jbeg: jend]
# Scalars
            self.psi_sep [jt] = data_block[jend]
            self.psi_axis[jt] = data_block[jend + 1]
            self.sgn_ip[jt] = int(data_block[jend + 2])
            self.Rmag    [jt] = data_block[jend + 3]
            self.Zmag    [jt] = data_block[jend + 4]

# Set units

        for attr in ('Rmesh', 'Zmesh', 'Rmag', 'Zmag'):
            self.__dict__[attr].units = '(m)'
        for attr in ('psi2d', 'psi1d', 'psi_sep', 'psi_axis'):
            self.__dict__[attr].units = '(Wb/rad)'
        for attr in ('rho2d', 'rho1d', 'q', 'sgn_ip'):
            self.__dict__[attr].units = ''
        self.vol.units  = '(m**3)'
        self.area.units = '(m**2)'
        self.F.units = '(T m)'

# Set description

        self.Rmesh.long_name  = 'R-mesh for equilibrium matrices'
        self.Zmesh.long_name  = 'Z-mesh for equilibrium matrices'
        self.Rmag.long_name   = 'R-coord of magnetic axis'
        self.Zmag.long_name   = 'Z-coord of magnetic axis'

        self.psi2d.long_name  = 'Psi matrix on cartesian R,z grid'
        self.rho2d.long_name  = 'rho matrix on cartesian R,z grid'

        self.psi1d.long_name  = 'Psi profile'
        self.rho1d.long_name  = 'rho profile'
        self.vol.long_name    = 'Volume within flux surfaces'
        self.area.long_name   = 'Area within flux surfaces'
        self.q.long_name      = 'Safety factor profile'
        self.F.long_name      = 'F function profile'

        self.sgn_ip.long_name   = 'Sign of plasma current'
        self.psi_axis.long_name = 'Psi at axis'
        self.psi_sep.long_name  = 'Psi at separatrix'

# Set dims

        self.Rmesh.dims = ('Rmesh', )
        self.Zmesh.dims = ('Zmesh', )
        self.psi2d.dims = ('time', 'Zmesh', 'Rmesh')
        self.rho2d.dims = ('time', 'Zmesh', 'Rmesh')
        for attr in ('Rmag', 'Zmag', 'psi_sep', 'psi_axis', 'sgn_ip'):
            self.__dict__[attr].dims = ('time',)
        for attr in ('psi1d', 'rho1d', 'vol', 'area', 'q', 'F'):
            self.__dict__[attr].dims = ('time', 'rho_eq')


class RABBIT_OUT():


    def __init__(self, nshot):

        """
        It reads all RABBIT binary output files for each NBI
        It generates a single NetCDF files concatenating all beams
        """

        nshot = int(nshot)
        self.shot = nshot
        rbpath = '%s/%d' %(rbhome, nshot)
        #rbpath = '/afs/ipp/home/d/dfajardo/rabbit/28053'
        beam = {}
        jnbi = 0
        fbin = '%s/beam%d/rtfi_result_oav.bin' %(rbpath, jnbi+1)
        if not(os.path.isfile(fbin)):
            print('Output file %s not found, exit' %fbin)
            return
        while(os.path.isfile(fbin)):
            beam[jnbi] = read_beam.read_outbeam(fbin)
            jnbi += 1
            fbin = '%s/beam%d/rtfi_result_oav.bin' %(rbpath, jnbi+1)
        nt, nv, nrho  = beam[0]['bdep'].shape
        nnbi = jnbi

        self.time    = beam[0]['time']
        self.rho_out = beam[0]['rho']
        self.dV      = beam[0]['dV']
        self.dArea   = beam[0]['dArea']

        self.time.dims    = ('time', )
        self.rho_out.dims = ('rho_out', )
        self.dV.dims    = ('time', 'rho_out')
        self.dArea.dims = ('time', 'rho_out')

        for sig, var in beam[0].items():
            if sig in ('time', 'rho', 'rabbit_version', 'dV', 'dArea'):
                continue
            darr = []
            for jnbi in range(nnbi):
                darr.append(beam[jnbi][sig])
            darr = np.array(darr, dtype=np.float32)
            self.__dict__[sig] = SFARR(np.sum(darr, axis=0))

            if darr.ndim == 2:
                self.__dict__[sig + 'Q'] = SFARR(np.transpose(darr, (1, 0)))
                self.__dict__[sig + 'Q'].dims = ('time', 'n_nbi')
                self.__dict__[sig].dims = ('time', )
            elif darr.ndim == 3:
                self.__dict__[sig + 'Q'] = SFARR(np.transpose(darr, (1, 2, 0)))
                self.__dict__[sig + 'Q'].dims = ('time', 'rho_out', 'n_nbi')
                self.__dict__[sig].dims = ('time', 'rho_out')
            elif darr.ndim == 4:
                self.__dict__[sig] = np.transpose(self.__dict__[sig], (0, 2, 1))
                self.__dict__[sig].dims = ('time', 'rho_out', 'spec')
                self.__dict__[sig + 'Q'] = SFARR(np.transpose(darr, (1, 3, 0, 2)))
                self.__dict__[sig + 'Q'].dims = ('time', 'rho_out', 'n_nbi', 'spec')

            self.__dict__[sig      ].units = var.units 
            self.__dict__[sig + 'Q'].units = var.units 
            self.__dict__[sig      ].long_name = var.long_name + ', sum_#NBI'
            self.__dict__[sig + 'Q'].long_name = var.long_name

        self.bdep_k1Q = np.sum(self.bdep_k1Q*self.bdepQ, axis=3)/np.sum(self.bdepQ, axis=3)
        self.bdep_k1Q.dims = ('time', 'rho_out', 'n_nbi')
        self.bdep_k1Q.units = beam[0]['bdep_k1'].units
        self.bdep_k1Q.long_name = beam[0]['bdep_k1'].long_name

# Neutron rate

        fneut = '%s/rtfi_nrate_oav.bin' %rbpath
        if os.path.isfile(fneut):
            self.nrate = read_beam.read_neut(fneut, nt, nrho)
            self.nrate.dims = ('time', 'rho_out')

# Derived signals

        self.Inbcd = np.sum(self.jnbcd*self.dArea, axis=1)
        self.Ifi   = np.sum(self.jfi  *self.dArea, axis=1)
        self.WfiPar  = np.sum(self.wfi_par *self.dV, axis=1)
        self.WfiParL = np.sum(self.wfi_parL*self.dV, axis=1)
        self.WfiPerp = np.sum(self.wfi_perp*self.dV, axis=1)
        #self.nyield  = np.sum(self.nrate  *self.dV, axis=1)

        self.Inbcd.units   = '(A)'
        self.Ifi.units     = '(A)'
        self.WfiPar.units  = '(J)'
        self.WfiParL.units = '(J)'
        self.WfiPerp.units = '(J)'
        #self.nyield.units  = '(cnt/s)'

        self.Inbcd.long_name   = 'NBI driven current'
        self.Ifi.long_name     = 'Fast ion current'
        self.WfiPar.long_name  = 'Parallel fast ion energy'
        self.WfiParL.long_name = 'Parallel fast ion energy, lab'
        self.WfiPerp.long_name = 'Perpendicular fast ion energy'
        #self.nyield.long_name  = 'Neutron rate'

        for sig in ('Inbcd', 'Ifi', 'WfiPar', 'WfiParL', 'WfiPerp'):#, 'nyield'):
            self.__dict__[sig].dims = ('time', )


if __name__ == '__main__':

    #import matplotlib.pylab as plt

    nshot = 33427
    #rb_out = RABBIT_OUT(nshot)
    rb_in = RABBIT_IN(nshot)
# Write shotfile
#    rb2sf(rb_in, rb_out)
# Write NetCDF
#    f_cdf = rb2cdf(rb_in, rb_out)
 
#    cv = netcdf.netcdf_file(f_cdf, 'r', mmap=False).variables
#    plt.plot(cv['time'].data, cv['pheatI'][:, 2], 'b-', label=r'P$_{heat, i}$')
#    plt.plot(cv['time'].data, cv['pheatE'][:, 2], 'r-', label=r'P$_{heat, e}$')
#    plt.legend()
    #wfi_par = rb_out.WfiPar
    #wfi_par_lab = rb_out.WfiParL
    #wfi_perp = rb_out.WfiPerp
    
    #print(wfi_par_lab)
    #print(wfi_perp)
    
    #print(wfi_par)
    #print(wfi_perp)
    #nrate = 
    #nt = len(rb_in.time)
    
    #print('difference between wfi_par and wfi_par_lab:')
    #print(wfi_par-wfi_par_lab)
    
    #wtor = rb_in.wtor
    
    #print('toroidal angular velocity')
    #print(wtor)
    
    
    print('ti==te?', rb_in.te == rb_in.ti)
    
    print('wtor', rb_in.wtor)

    #plt.subplot(1, 2, 1)
    #plt.plot(rb_in.rho_in, rb_in.te[nt//2, :])
    #plt.subplot(1, 2, 2)
    #plt.plot(rb_in.rho1d[nt//2, :], rb_in.vol[nt//2, :])

    #plt.show()
