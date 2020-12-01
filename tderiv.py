import numpy as np
from scipy.optimize import minimize


def func(x, *p):
    a, b, c = p
    return a - b*c*np.exp(-x/c)


def exp_fit(xin, yin, lder):
    p_init = [ yin[0], lder, 0.1]
    bounds = [(0., 1.e3), (lder, lder), (1e-6, 10.)]
    err = lambda p: np.mean((func(xin, *p) - yin)**2)
    p_opt = minimize( \
        err, \
        p_init, \
        bounds=bounds, \
        options={'maxiter':100000}, \
        method="L-BFGS-B").x
    return p_opt


def half_der(time, val, jt1, jt2):

    jt2 = min(jt2, len(val)-1)
    tfit = time[jt1:jt2+1]
    yfit = val[jt1:jt2+1]
    n_dim = jt2 - jt1 + 1
    if n_dim <= 1:
        return 0
    elif n_dim == 2:
        return (val[jt2]-val[jt1])/(time[jt2]-time[jt1])
    else:
        dtj = time[jt1: jt2+1] - time[jt1]
        dvj = val[jt1: jt2+1] - val[jt1]
        sx  = np.sum(dtj)
        sy  = np.sum(dvj)
        sxx = np.sum(dtj**2)
        sxy = np.sum(dtj*dvj)
        return (n_dim*sxy - sx*sy)/(n_dim*sxx - sx**2)


def dvdt_ic(val, ref, dt=0.001, dtint=0.0301, dt_ic=0.0301, err_ref=1.e5, err_abs=1e8):

# git solve dW/dt smearring issue?
#    dtint = 0.005
#    dt_ic = 0.005

    n_tim = len(val)
    time = dt*np.arange(1, n_tim+1)
    dref = np.abs(np.diff(ref))
    ndt = int(dtint/dt)
    print('tderv:dvdt_ic', n_tim, time[-1])

    istat = np.zeros(n_tim, dtype=np.int32)
    li = np.zeros(n_tim+ndt, dtype=np.int32)
    ri = np.zeros(n_tim+ndt, dtype=np.int32)
    dvdt  = np.zeros(n_tim, dtype = np.float32)
    dvdt2 = np.zeros(n_tim, dtype = np.float32)
    kl1 = 0

    if np.abs(ref[1]-ref[0]) > err_ref :
        istat[0] = -1
    else:
        istat[0] = 1
    if np.abs(ref[-1]-ref[-2]) > err_ref :
        istat[-1] = -1
    else:
        istat[-1] = 2
    for jt in range(1, n_tim-1):
        if dref[jt] > err_ref:
            if dref[jt-1] > err_ref:
                istat[jt] = -1
            else:
                istat[jt] = 2
        elif dref[jt-1] > err_ref:
            istat[jt] = 1
        else:
            istat[jt] = 0

    int1flg = False

    for jt in range(n_tim):
        if istat[jt] == -1:
            li[jt] = 1
            ri[jt] = 1
        elif istat[jt] == 1:
            for j in range(jt+2, n_tim+1):
                li[j-1] = j - jt - 1
                if istat[j-1] == 2:
                    break
            jr = j
            for j in range(jt+1, jr+1):
                ri[j-1] = jr - j
    
            jright = 0
            if int1flg:
                tfit = np.zeros(jr - jt)
                vfit = np.zeros(jr - jt)
                for j in range(jt+1, jr+1):
                    t_loc = time[j-1] - time[jt]
                    tfit[j-jt-1] = t_loc
                    vfit[j-jt-1] = val[j-1]
                    if t_loc > dtint:
                        j += -1
                        break
                jright = j+1
                n_fit = jright - jt - 1
                if kl1 > 0:
                    kl1 += -1
                    il = kl1 - li[kl1-1]
                    l_der = half_der(time, val, il-1, kl1-1)
                else:
                    kl1 *= -1
                if jt > 0:
                    if jright > jt+2:
                        lder = 1e-6*(l_der + ref[jt] - ref[kl1-1])
                        par = exp_fit(tfit[:n_fit], 1e-6*vfit[:n_fit], lder)
                        dvdt[jt:jright-1] = 1e6*lder*np.exp(-tfit[:n_fit]/par[2])
                else:
                    dvdt[jt:jright-1] = 0.

                istat[jt:jt+n_fit] += 3 # Exclude discotinuities (+ interval) from integral correction

                dvdt[kl1-1] = l_der
                dvdt[kl1:jt] = l_der + ref[kl1:jt] - ref[kl1-1]
                istat[kl1:jt] = 6

            for i_test in range(jright, jr+1):
                if (li[i_test-1] > 1) and (istat[i_test-1] < 3):
                    for j in range(2, li[i_test-1]+1):
                        if time[i_test-1]-time[i_test-j-1] > dtint: 
                            j += -1
                            break
                    li[i_test-1] = min([j, li[i_test-1] - 1])
                if (ri[i_test-1] > 1) and (istat[i_test-1] < 3):
                    for j in range(2, ri[i_test]+1):
                        if time[i_test+j-1] - time[i_test-1] > dtint:
                            j += -1
                            break
                    ri[i_test-1] = min([j, ri[i_test-1] - 1])

            jt = jr - 1
            kl3 = abs(kl1)
            kl1 = jr
            if jright > jr:
                kl1 = -jright + 1
                l_der = dvdt[jright-2]
            if (not int1flg):
                kl2 = abs(kl1)
            int1flg = True

    li[0]  = 0
    ri[-1] = 0

    for jt in range(n_tim):
        if istat[jt] < 3:
            l_der = half_der(time, val, jt - li[jt], jt)
            r_der = half_der(time, val, jt, jt + ri[jt])
            if np.abs(r_der - l_der) < err_abs:
                dvdt[jt] = (ri[jt]*r_der + li[jt]*l_der)/(ri[jt] + li[jt])

# Integral correction

    ndt_ic = int(dt_ic/dt)
    int_dv = np.zeros(n_tim, dtype=np.float32)
    dvn    = np.zeros(n_tim, dtype=np.float32)

    if (kl3 > kl2) and (kl2 > 1) and (kl3 < n_tim):
        int_dv[kl2-1] = val[kl2-1]
        imax = min(kl3+2*ndt_ic, n_tim-1)
        ind = range(kl2, imax)
        for jt in ind:
# Time integral of dvdt
            int_dv[jt] = int_dv[jt-1] + (dvdt[jt-1] + dvdt[jt])*dt*0.5

        dvn[ind] = val[ind] - int_dv[ind]
   
        i0 = max(kl2-ndt_ic, 0)
        i1 = min(kl3+2*ndt_ic+6, n_tim)
        for jt in range(i0, i1-1):
            il = max(jt+1-ndt_ic, i0+1)
            ir = min(jt+1+ndt_ic, i1-1)
            int_dv[jt] = np.average(dvn[il:ir+1])

# Calculate d/dt of the time integral
        print('GIT test', i0, i1-2, kl3+2*ndt_ic+6)
        ind = range(i0, i1-2)
        dvdt2[ind] = dvdt[ind] + np.gradient(int_dv[ind], dt)
    else:
        print('no integral correction in routine time_der4')

    return dvdt2


def corr_u(time, flux_in, RC):

# Time derivative with i*omega/(1+i*RC*omega) filter (low pass)

    nt = len(flux_in)
    dt = time[1] - time[0]
    omega = 2*np.pi*np.fft.fftfreq(nt, d=dt)
    print(time[-1] - time[0])
    print('omega=', omega)
    print('Fourier transform of Phi(B)')
    flux_fft = np.fft.fft(flux_in)
    fac = 1.j*omega/(1 + 1.j*omega*RC) # d/dt with low pass filter
                                        # d/dt only for RC=0
    tmp = fac*flux_fft
    print('Inverse Fourier transform, RC -> low pass filter')
    flux_out  = np.fft.ifft(tmp)
    rflux_out = np.real(flux_out)
#    rflux_out[:10] = 0.

    return rflux_out
