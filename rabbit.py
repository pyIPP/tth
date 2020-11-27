#!/usr/bin/env python

import os, sys
sys.path.append('/afs/ipp/home/g/git/trview')
sys.path.append('/afs/ipp/home/g/git/python/rbview')
import tr_batch
import rb_io, write_rabbit


if os.getenv('TRVIEW_RBHOME') is None:
    rbhome = '/toks/work/%s/rabbit' %os.getenv('USER')
else:
    rbhome = os.getenv('TRVIEW_RBHOME')

#RABBIT executable path
rb_exec_path = '/afs/ipp-garching.mpg.de/home/m/markusw/fortran/IDAbesrad/branches/d3dcopy3.3/build/rabbit'


class RABBIT():


    def __init__(self, nshot, tbeg, tend):

        '''
        creates RABBIT input for shotfile number nshot
        runs RABBIT
        reads the RABBIT output, in particular wfi_par, wfi_per
        '''

        rbpath = '/toks/work/%s/rabbit/%d' %(os.getenv('USER'), nshot)
        #rbpath = '/afs/ipp/home/d/dfajardo/rabbit/28053/'
        os.system('mkdir -p %s' %rbpath)
    
#------------------------get experimental data-----------------------------
# Get settings
        gpar = tr_batch.GPAR(nshot, tbeg, tend)
# Get data
        self.rb_in = tr_batch.INPUT(gpar, code='rb')
# Write RABBIT input
        write_rabbit.write_rabbit(self.rb_in.profiles, self.rb_in.equ, self.rb_in.sig.zef, self.rb_in.nbi, self.rb_in.ion_d)

#---------------------------execute RABBIT---------------------------------
        print('Executing RABBIT...')
        os.system('%s %s' %(rb_exec_path, rbpath))
        print('Finished executing RABBIT!')

#----------------------read RABBIT output files----------------------------    
        self.rb_out = rb_io.RABBIT_OUT(nshot)


if __name__ == '__main__':


    #import numpy as np
    #import matplotlib.pylab as plt
    #import argparse

    #parser = argparse.ArgumentParser(description='RABBIT stand-alone')
    #parser.add_argument('-s', '--shot', type=int, help='Shot number' , required=False)
    #parser.add_argument('-b', '--tbeg', type=float, help='Initial time', required=False)
    #parser.add_argument('-e', '--tend', type=float, help='End time', required=False)
   
    #args = parser.parse_args()

    #if args.shot is None:
    #    shot = 28053
    #else:
    #    shot = args.shot
    #if args.tbeg is None:
    #    tbeg = 0.
    #else:
    #    tbeg = args.tbeg
    #if args.tend is None:
    #    tend = 10.
    #else:
    #    tend = args.tend
    shot = int(input('nshot: ')) #33500
    tbeg = 0.#0.001
    tend = 10.#5.511
    rb = RABBIT(shot, tbeg, tend)

    #wfi_par  = np.sum(rb.rb_out.wfipar * rb.rb_out.dV, axis = 1)
    #wfi_perp = np.sum(rb.rb_out.wfiperp* rb.rb_out.dV, axis = 1)
    wfi_par = rb.rb_out.WfiPar
    wfi_par_lab = rb.rb_out.WfiParL
    wfi_perp = rb.rb_out.WfiPerp
    
    print('Wfi par',wfi_par, 'shape = ', wfi_par.shape)
    print('Wfi perp', wfi_perp, 'shape = ', wfi_perp.shape)
    
    c_wfi_par = 1.5
    c_wfi_perp = 0.75
    wfi = c_wfi_par*wfi_par_lab + c_wfi_perp*wfi_perp
    
    print('Wfi = ', wfi, 'shape = ', wfi.shape)
    
    
    nt = len(rb.rb_out.time)
    print('nt = ', nt)#, 'Should be 541')
    
    
    print('difference between wfi_par and wfi_par_lab:')
    print(wfi_par-wfi_par_lab)
    
    wtor1 = rb.rb_in.profiles.prof['V'].fit_rhop
    print('toroidal angular velocity')
    print(wtor1)

    wtor = rb.rb_in.wtor
    
    print('toroidal angular velocity')
    print(wtor)

    #plt.plot(rb.rb_out.time, wfi_par , 'b-', label='Par')
    #plt.plot(rb.rb_out.time, wfi_perp, 'r-', label='Perp')
    #plt.legend()
    #plt.show()
