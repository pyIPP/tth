"""
    Module for the calculation of the NBI power losses
    (function for single NBI calculation)
"""
__author__='Marco Cavedon'
__date__='14 January 2014'
__version__='2.1'

import numpy as np

# Tabulated values of beam stopping sigma

EdivA_tab = [5.0, 8.0, 10., 15., 20., 30., 40.,\
             50., 55., 60., 80., 100.]
sig_tab = [1.490, 1.290, 1.200, 1.040, 0.895, 0.640,\
           0.458, 0.352, 0.318, 0.292, 0.230, 0.196]


def riploss(nshot, jnbi):

# Rip Losses/P_NBI

    if (nshot < 14049):
        pars = 2*[0.07, 0.0, 0.0, 0.07]
    else:
        pars = [0.07, 0.0, 0.0, 0.07, \
                0.00, 0.0, 0.0, 0.00]
    return pars[jnbi]


def solloss(nshot, fraci, EdivA, jnbi):

    """ SOL Losses/P_NBI/ne_sol
    input: - EdivA [keV/AMU]: beam energy over
             mass of the injected particles
           - fraci [%]: power distribution of the NBI
           - jnbi: NBI #source (INJ1:0-3, INJ2:4-7)
    """
    
    if (nshot < 14049):
        pars  = 2*[0.5995, 0.6463, 0.6038, 0.5642]
    else:
        pars  = [1.019, 1.099, 1.026, 0.959, \
                 1.099, 1.700, 1.360, 1.026]

    solloss = 0.
    for jmix in range(len(fraci)):
        sigma = np.interp(EdivA/(jmix+1), EdivA_tab, sig_tab)
        solloss += fraci[jmix]*sigma

    return pars[jnbi]*0.06*solloss


class FAFPAR:


    """
    NBI losses following parametriation
    made on data from FAFNER (McCarthy, Staebler).
    For the definitions of the losses see:
    http://www.ipp.mpg.de/~git/tot/index.html
    """

    def __init__(self, nshot, fraci, ne_lav, Te_lav, ne_sep, tau_sd, EdivA, jnbi):

        # NBI losses parameters
        pdata1 = {}
        pdata2 = {}
        pdata1['wfi']      = 2*[922.5, 972.7, 1006.2, 944.8]
        pdata1['shine']    = 2*[1.49, 1.16, 1.10, 1.42]
        pdata1['cxloss1']  = 2*[2.836, 5.505, 2.873, 2.344]
        pdata1['cxloss2']  = 2*[1.189, 1.414, 1.262, 1.130]
        pdata1['orbloss1'] = 2*[0.0464, 0.0395, 0.0349, 0.0432]
        pdata1['orbloss2'] = 2*[0.492, 0.417, 0.496, 0.514]
        pdata1['orbloss3'] = 8*[0.0]

        pdata2['wfi'] = pdata1['wfi']
        pdata2['shine']   = [0.894, 0.696, 0.66, 0.852, \
                             0.84 , 0.54 , 0.48, 0.78]
        pdata2['cxloss1'] = 8*[0.5]
        pdata2['cxloss2'] = 8*[1.2]
        pdata2['orbloss1'] = [0.0464, 0.0395, 0.0349, 0.0432, \
                              0.0395, 0.0454, 0.015 , 0.0349]
        pdata2['orbloss2'] = [0.45, 0.45, 0.45, 0.45, \
                              0.0 , 0.0 , 0.0 , 0.0]
        pdata2['orbloss3'] = [0.0, 0.0, 0.0, 0.0, \
                              0.0, 0.7, 1.3, 0.0]

        if (nshot < 14049):
            pdata = pdata1
        else:
            pdata = pdata2

# Shine-Through/P_NBI

        self.st = np.zeros_like(ne_lav)
        for jmix in range(len(fraci)):
            sigma = np.interp(EdivA/(jmix+1), EdivA_tab, sig_tab)
            self.st += fraci[jmix]*np.exp(-ne_lav*sigma/pdata['shine'][jnbi])

# Orbit losses/P_NBI

        self.orb = pdata['orbloss1'][jnbi]* \
                   np.power(ne_sep*Te_lav, pdata['orbloss2'][jnbi]) + \
                   pdata['orbloss3'][jnbi]*self.st

# CS losses/P_NBI

        self.cx = pdata['cxloss1'][jnbi]* \
                  np.power(tau_sd, pdata['cxloss2'][jnbi])
    
# Fast ion energy [J]
        self.wfi = pdata['wfi'][jnbi]*tau_sd*4.e-4



class TRAPAR:

    
    def __init__(self, nshot, ipl_sign, Te_lav, ne_lav, Ip, EdivA, pf, Zeff, jnbi):

        """
        NBI losses following a new parametrization 
        made on data from TRANSP/NUBEAM (Cavedon, Tardini, 2012)
        ipl_sign='straight' -> Conventional Magnetic Configuration
        ipl_sign='rev' -> Reversed Magnetic Configuration
        """

        pdata = {}

        for channel in ('ST', 'ORB', 'CXint', 'CXext', 'WFI'):
            pdata[channel] = {}
            for field in ('straight', 'rev'):
                pdata[channel][field] = {}

        pdata['ST']['straight']['ne'] = [1.019, 0.9832, 0.9782, 1.0252, \
                                         1.003, 0.7732, 0.7616, 1.006]
        pdata['ST']['straight']['E'] = [-0.5302, -0.513 , -0.5078, -0.53, \
                                        -0.5231, -0.4161, -0.4113, -0.521]
        pdata['ST']['straight']['Te'] = [-0.0765, -0.0743, -0.0723, -0.0754, \
                                         -0.076 , -0.0567, -0.0559, -0.073]
        pdata['ST']['straight']['Zeff'] = [0.1187, 0.1155, 0.1163, 0.1215, \
                                           0.1139, 0.0883, 0.0907, 0.117]
        pdata['ST']['straight']['pf'] = [-0.1414, -0.1632, -0.0517, -0.0376, \
                                         -0.2222, -0.3629, -0.2806, -0.068]
        pdata['ST']['straight']['C'] = [0.2626, 0.2191, 0.2208, 0.2674, \
                                        0.2329, 0.1373, 0.1182, 0.235]

        pdata['ST']['rev']['ne'] = [1.0163, 0.9788, 0.9737, 1.0219, \
                                    0.9995, 0.7703, 0.7622, 1.004]
        pdata['ST']['rev']['E']  = [-0.5304, -0.5147, -0.5086, -0.5302, \
                                      -0.5243, -0.4185, -0.4114, -0.522]
        pdata['ST']['rev']['Te'] = [-0.0761, -0.0732, -0.0713, -0.0752, \
                                    -0.0752, -0.0555, -0.0559, -0.073]
        pdata['ST']['rev']['Zeff'] = [0.1231, 0.1216, 0.1244, 0.1268, \
                                      0.1189, 0.0907, 0.0959, 0.124]
        pdata['ST']['rev']['pf'] = [-0.1431, -0.1639, -0.0543, -0.0405, \
                                    -0.2231, -0.3596, -0.2802, -0.071]
        pdata['ST']['rev']['C'] = [0.2619, 0.2171, 0.2198, 0.2665, \
                                   0.2313, 0.1359, 0.1188, 0.234]

        pdata['ORB']['straight']['C'] = [0.0009, 0.0005, 0.0005, 0.0008, \
                                         0.0005, 0.0001, 0.0001, 0.000]
        pdata['ORB']['straight']['Te'] = [0.6456, 0.6985, 0.6985, 0.6389, \
                                          0.6923, 0.727 , 0.7948, 0.688]
        pdata['ORB']['straight']['ne'] = [-0.2042, -0.7207, -0.7002, -0.1458, \
                                          -0.4589, -1.3537, -1.456 , -0.426]
        pdata['ORB']['straight']['Ip'] = [-2.5337, -2.6013, -2.6177, -2.4989, \
                                          -2.6486, -3.03  , -3.1295, -2.633]
        pdata['ORB']['straight']['E'] = [0.763 , 1.0394, 1.0352, 0.7613, \
                                         0.9274, 1.5761, 1.5166, 0.922]
        pdata['ORB']['straight']['pf'] = [-0.1852, 0.052, -0.1187, -0.2942, \
                                           0.0806, 0.9632, 0.5579, -0.154]
        pdata['ORB']['straight']['Zeff'] = [0.6165, 0.5408, 0.5458, 0.6254, \
                                            0.5877, 0.4235, 0.4969, 0.598]

        pdata['ORB']['rev']['C'] = [0.0088, 0.0367, 0.0422, 0.0113, \
                                    0.0272, 0.1025, 0.0748, 0.033]
        pdata['ORB']['rev']['Te'] = [0.0804, 0.0023, -0.0177, 0.0544, \
                                     0.0258, 0.0273,  0.0043, -0.011]
        pdata['ORB']['rev']['ne'] = [0.7886, 0.816, 0.821 , 0.8122, \
                                     0.8133, 0.714, 0.7868, 0.824]
        pdata['ORB']['rev']['Ip'] = [-1.6839, -1.2256, -1.094 , -1.5, \
                                     -1.3716, -0.6652, -0.8154, -1.162]
        pdata['ORB']['rev']['E'] = [0.3929, 0.1749,  0.1351, 0.3259, \
                                    0.2311, 0.0653, -0.0021, 0.167]
        pdata['ORB']['rev']['pf'] = [-0.3812, -0.4589, -0.4781, -0.4177, \
                                     -0.4349, -0.3585, -0.4279, -0.469]
        pdata['ORB']['rev']['Zeff'] = [0.264 , 0.1627, 0.1533, 0.2499, \
                                       0.1765, 0.1335, 0.1681, 0.165]

        pdata['CXint']['straight']['C'] = [1.5788, 1.5608, 1.7595, 1.5814, \
                                           1.2892, 1.0385, 1.1255, 1.464]
        pdata['CXint']['straight']['E'] = [-0.7779, -0.7115, -0.6601, -0.7201, \
                                           -0.7153, -0.616 , -0.5696, -0.655]
        pdata['CXint']['straight']['Zeff'] = [-0.9361, -0.8599, -0.839 , -0.9119, \
                                              -0.8867, -0.8158, -0.8152, -0.863]
        pdata['CXint']['straight']['ne'] = [-1.5631, -1.9837, -2.148 , -1.6879, \
                                            -1.7945, -2.1479, -2.3056, -1.992]
        pdata['CXint']['straight']['Te'] = [0.8233, 0.9433, 0.9146, 0.8086, \
                                            0.9201, 1.1364, 1.1021, 0.896]
        pdata['CXint']['straight']['Ip'] = [0.1391, 0.1109,  0.0576, 0.0938, \
                                            0.1378, 0.0216, -0.0282, 0.07]
        pdata['CXint']['straight']['pf'] = [0.0883, 0.2169, 0.1418, 0.0571, \
                                            0.2334, 0.6999, 0.5427, 0.150]

        pdata['CXint']['rev']['C'] = [1.702 , 2.1279, 2.1862, 1.8438, \
                                      1.6458, 0.9258, 0.8109, 1.951]
        pdata['CXint']['rev']['E'] = [-0.8771, -0.8934, -0.8122, -0.8432, \
                                      -0.8673, -0.6761, -0.6429, -0.816]
        pdata['CXint']['rev']['Zeff'] = [-0.9371, -0.8798, -0.8793, -0.9366, \
                                         -0.8963, -0.7362, -0.7567, -0.885]
        pdata['CXint']['rev']['ne'] = [-1.4229, -1.8436, -2.0462, -1.5812, \
                                       -1.6595, -1.9613, -1.895 , -1.91]
        pdata['CXint']['rev']['Te'] = [0.8705, 0.9839, 0.9546, 0.865, \
                                       0.9609, 1.0909, 1.0217, 0.938]
        pdata['CXint']['rev']['Ip'] = [0.6079, 0.8143, 0.7093, 0.5728, \
                                       0.8311, 0.6008, 0.5595, 0.689]
        pdata['CXint']['rev']['pf'] = [0.0441, 0.1887, 0.1111, -0.0067, \
                                       0.1928, 0.3601, 0.0642,  0.056]

        pdata['CXext']['straight']['C'] = [0.0268, 0.0316, 0.031, 0.028, \
                                           0.0268, 0.025, 0.0138, 0.026]
        pdata['CXext']['straight']['E'] = [-0.4023, -0.4553, -0.4407, -0.3999, \
                                           -0.4310, -0.0338, -0.1049, -0.4120]
        pdata['CXext']['straight']['Zeff'] = [0.2995, 0.3435, 0.2885, 0.2585, \
                                              0.3440, 0.2918, 0.2787, 0.2800]
        pdata['CXext']['straight']['ne'] = [-0.0276, -0.4666, -0.4560, -0.0346, \
                                            -0.3014, -1.4345, -1.1162, -0.298]
        pdata['CXext']['straight']['Te'] = [0.9534, 1.1309, 1.0857, 0.9088, \
                                            1.0889, 1.2482, 1.2557, 1.0220]
        pdata['CXext']['straight']['Ip'] = [-0.4423, -0.3888, -0.3536, -0.3864, \
                                            -0.4527, -0.5940, -0.6405, -0.3700]
        pdata['CXext']['straight']['pf'] = [-0.0322, 0.0782, 0.0310, -0.0456, \
                                             0.0879, 0.6933, 0.4494,  0.0370]

        pdata['CXext']['rev']['C'] = [0.2036, 0.2920, 0.3066, 0.2114, \
                                      0.2627, 0.3907, 0.3065, 0.280]
        pdata['CXext']['rev']['E'] = [-0.7219, -0.9119, -0.9794, -0.768, \
                                      -0.8321, -0.8859, -0.9236, -0.920]
        pdata['CXext']['rev']['Zeff'] = [0.0500, 0.0934, 0.1103, 0.0556, \
                                         0.0748, 0.1607, 0.1998, 0.092]
        pdata['CXext']['rev']['ne'] = [-0.0122, -0.1795, -0.1556, -0.0287, \
                                       -0.1709, -0.2764, -0.0623, -0.148]
        pdata['CXext']['rev']['Te'] = [0.7448, 0.7628, 0.7454, 0.7422, \
                                       0.7607, 0.5999, 0.5216, 0.7400]
        pdata['CXext']['rev']['Ip'] = [0.3757,  0.5364, 0.4713, 0.3637, \
                                       0.5282, -0.1206, 0.1430, 0.4770]
        pdata['CXext']['rev']['pf'] = [-0.3273, -0.1596, -0.2316, -0.370, \
                                       -0.1748, -0.0464, -0.2821, -0.269]

        pdata['WFI']['straight']['C'] = [4333.4475, 4120.8556, 4259.2208, 4421.2517, \
                                         3997.6816, 4731.2498, 4500.7905, 4073.35]
        pdata['WFI']['straight']['ne'] = [-0.8676, -0.9229, -0.9662, -0.899, \
                                          -0.8809, -0.9441, -1.0554, -0.929]
        pdata['WFI']['straight']['E'] = [0.6266, 0.6873, 0.7134, 0.6528, \
                                         0.6535, 0.5614, 0.6786, 0.6930]
        pdata['WFI']['straight']['Te'] = [0.7204, 0.7235, 0.7082, 0.7029, \
                                          0.7312, 0.8516, 0.8032, 0.7090]
        pdata['WFI']['straight']['Ip'] = [0.0345, -0.0137, -0.0098,  0.0353, \
                                         -0.0084,  0.0286, -0.0185, -0.0120]
        pdata['WFI']['straight']['pf'] = [0.0800, 0.0426, 0.0265, 0.0599, \
                                          0.0721, 0.3278, 0.2451, 0.0500]

        pdata['WFI']['rev']['C'] = [8118.7247, 7643.8437, 6013.1903, 6914.7088, \
                                    8638.8638, 6995.6925, 4252.9034, 6043.776]
        pdata['WFI']['rev']['ne'] = [-0.9941, -1.1924, -1.1995, -1.0214, \
                                     -1.1338, -1.4925, -1.2878, -1.146]
        pdata['WFI']['rev']['E'] = [0.3946, 0.4653, 0.5851, 0.4859, \
                                    0.3809, 0.4512, 0.6762, 0.5520]
        pdata['WFI']['rev']['Te'] = [0.7759, 0.8167, 0.7874, 0.7557, \
                                     0.8103, 0.9124, 0.8874, 0.783]
        pdata['WFI']['rev']['Ip'] = [0.4476, 0.4753, 0.3301, 0.3521, \
                                     0.5687, 0.4443, 0.1326, 0.3410]
        pdata['WFI']['rev']['pf'] = [0.4302, 0.4583, 0.3400, 0.3578, \
                                     0.5073, 0.7273, 0.4127, 0.3570]

        if (nshot < 14049):
            print('No TRANSP calculations for old NBI geometry,  i.e. #shot < 14049')
            return

        ind_pos = (ne_lav > 0) & (Te_lav > 0) & (Ip > 0) & (pf > 0)

# Shine-Through/P_NBI

        pars = pdata['ST'][ipl_sign]
        self.st = np.zeros_like(ne_lav)
        self.st[ind_pos] = np.exp( -1./pars['C'][jnbi]* \
            np.power(ne_lav[ind_pos],  pars['ne'  ][jnbi])* \
            np.power(EdivA,            pars['E'   ][jnbi])* \
            np.power(Te_lav[ind_pos],  pars['Te'  ][jnbi])* \
            np.power(pf[ind_pos],      pars['pf'  ][jnbi])* \
            np.power(Zeff[ind_pos],    pars['Zeff'][jnbi]) )

# Orbit, CXext, CXint/P_NBI

        self.cx_d = {}
        for lbl in ('ORB', 'CXext', 'CXint'):
            pars = pdata[lbl][ipl_sign]
            self.cx_d[lbl] = np.zeros_like(ne_lav)
            self.cx_d[lbl][ind_pos] = pars['C'][jnbi]*( \
                np.power(ne_lav[ind_pos], pars['ne'][jnbi])* \
                np.power(Te_lav[ind_pos], pars['Te'][jnbi])* \
                np.power(Ip[ind_pos],     pars['Ip'][jnbi])* \
                np.power(EdivA,           pars['E'][jnbi])* \
                np.power(pf[ind_pos],     pars['pf'][jnbi])* \
                np.power(Zeff[ind_pos],   pars['Zeff'][jnbi]) )

        self.cx = self.cx_d['CXint'] + self.cx_d['CXext']
        self.orb = self.cx_d['ORB']

# Fast Ion Energy [J]
# Unit is J despite 1e-6 pre-factor

        pars = pdata['WFI'][ipl_sign]
        self.wfi = np.zeros_like(ne_lav)
        self.wfi[ind_pos] = 1e-6*pars['C'][jnbi]* \
            np.power(EdivA,           pars['E'][jnbi])* \
            np.power(ne_lav[ind_pos], pars['ne'][jnbi])* \
            np.power(Te_lav[ind_pos], pars['Te'][jnbi])* \
            np.power(Ip[ind_pos],     pars['Ip'][jnbi])* \
            np.power(pf[ind_pos],     pars['pf'][jnbi])
