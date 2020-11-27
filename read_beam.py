import os
from bin2arr import bin2sfarr


def read_outbeam(filename):


    struct_len = 4 #struct.calcsize('f')

    beam = {}
    psigs = ('pheatI', 'pheatE', 'pheat', 'pshine', 'prot', 'ploss', 'pcx')

    with open(filename, mode='rb') as f:

        file_size = os.path.getsize(filename)

        ntime = int(bin2sfarr(f))
        nrho  = int(bin2sfarr(f))
        nv    = int(bin2sfarr(f))

        nlen = nrho*ntime*struct_len

        beam['time'] = bin2sfarr(f, (ntime, ))
        beam['rho' ] = bin2sfarr(f, (nrho, ))

# The sequence is important!
        for prof2 in ('bdens', 'press', 'powe', 'powi', 'jfi', 'jnbcd', 'dV'):
            beam[prof2] = bin2sfarr(f, (ntime, nrho))

        for prof3 in ('bdep', 'bdep_k1'):
            beam[prof3] = bin2sfarr(f, (ntime, nv, nrho))

        for sig in psigs:
            beam[sig] = bin2sfarr(f, (ntime, ))

# test to see if torque calculations have been preformed
        if file_size - f.tell() > nlen:
            rabbit_version_strlen = bin2sfarr(f, fmt='i')
            if int(rabbit_version_strlen) > 0:
                beam['rabbit_version'] = bin2sfarr(f, (rabbit_version_strlen, ), fmt='s')
#                print(beam['rabbit_version'])
                beam['dArea'  ] = bin2sfarr(f, (ntime, nrho))
                beam['torqdepo'] = bin2sfarr(f, (ntime, nv, nrho))
                beam['torqjxb'] = bin2sfarr(f, (ntime, nv, nrho))
                beam['torqe'  ] = bin2sfarr(f, (ntime, nrho))
                beam['torqi'  ] = bin2sfarr(f, (ntime, nrho))
            else:
                beam['rabbit_version'] = 'undef'

        next_line = f.read(struct_len) # nleg
        if not next_line:
            return beam

# test to see if wfi calculations have been preformed
        if file_size - f.tell() >= nlen:
            beam['wfi_par' ] = bin2sfarr(f, (ntime, nrho))
            beam['wfi_perp'] = 1.5*beam['press'] - beam['wfi_par']

# test to see if wfi_par_lab calculations have been preformed
        if file_size - f.tell() >= nlen:
            beam['wfi_parL'] = bin2sfarr(f, (ntime, nrho))

# Variables' description

    beam['time'    ].long_name = 'Time'
    beam['rho'     ].long_name = 'rho'
    beam['bdens'   ].long_name = 'Fast-ion Density'
    beam['press'   ].long_name = 'Fast-ion Pressure'
    beam['powe'    ].long_name = 'Power Density Profile to Electrons'
    beam['powi'    ].long_name = 'Power Density Profile to Ions'
    beam['jfi'     ].long_name = 'Fast-ion Current Density'
    beam['jnbcd'   ].long_name = 'Driven Current Density'
    beam['bdep'    ].long_name = 'Particle Source Density (per Energy Component)'
    beam['bdep_k1' ].long_name = r'Particle Source Average Pitch $\langle v_\parallel/v\rangle$ (per Energy Component)'
    beam['dV'      ].long_name = 'Volume of Radial Cell'
    beam['pheatI'  ].long_name = 'Ion Heating Power'
    beam['pheatE'  ].long_name = 'Electron Heating Power'
    beam['pheat'   ].long_name = 'Heating Power'
    beam['pshine'  ].long_name = 'Shine-thru Power'
    beam['prot'    ].long_name = 'Power to Rotation'
    beam['ploss'   ].long_name = 'Beam Power Losses in SOL'
    beam['pcx'     ].long_name = 'Charge Exchange Loss'
    beam['dArea'   ].long_name = 'Poloidal Cross-Sectional Area of Radial Cell'
    beam['torqdepo'].long_name = 'Deposited Torque Density'
    beam['torqjxb' ].long_name = 'JxB Torque Density'
    beam['torqe'   ].long_name = 'Collisional Torque Density to Electrons'
    beam['torqi'   ].long_name = 'Collisional Torque Density to Ions'
    beam['wfi_par' ].long_name = 'Stored Fast-Ion Energy Density (parallel, plasma frame)'
    beam['wfi_perp'].long_name = 'Stored Fast-Ion Energy Density (perpendicular)'
    beam['wfi_parL'].long_name = 'Stored Fast-Ion Energy Density (parallel, lab frame)'

    beam['time'    ].units = r'(s)'
    beam['rho'     ].units = r'(-)'
    beam['bdens'   ].units = r'(# m$^{-3}$)'
    beam['press'   ].units = r'(Pa)'
    beam['powe'    ].units = r'(W / m$^3$)'
    beam['powi'    ].units = r'(W / m$^3$)'
    beam['jfi'     ].units = r'(A / m$^2$)'
    beam['jnbcd'   ].units = r'(A / m$^2$)'
    beam['bdep'    ].units = r'(# m$^{-3}$ s$^{-1}$)'
    beam['bdep_k1' ].units = r'(-)'
    beam['dV'      ].units = r'(m$^3$)'
    for sig in psigs:
        beam[sig].units = r'(W)'
    beam['dArea'   ].units = r'(m$^2$)'
    beam['torqdepo'].units = r'(Nm / m$^3$)'
    beam['torqjxb' ].units = r'(Nm / m$^3$)'
    beam['torqe'   ].units = r'(Nm / m$^3$)'
    beam['torqi'   ].units = r'(Nm / m$^3$)'
    beam['wfi_par' ].units = r'(J / m$^3$)'
    beam['wfi_perp'].units = r'(J / m$^3$)'
    beam['wfi_parL'].units = r'(J / m$^3$)'

    return beam


def read_neut(filename, ntime, nrho):

    with open(filename, mode='rb') as f:
       neut = bin2sfarr(f, (ntime, nrho))
       neut.long_name = 'Neutron rate per unit volume'
       neut.units =  r'(# m$^{-3}$ s$^{-1}$)'
    return neut
