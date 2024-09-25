# import libraries
import numpy as np
import utilities as ut
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import random


def siDHM_forloop(holo1, holo2, wavelength, dxy):
    # Function to retrieve phase maps of off-axis SI-DHM via for loop.
    # inputs:
    # holo1 - [1 of 2] phase-shift holograms
    # holo2 - [2 of 2] phase-shift holograms
    # wavelength - Wavelength of the illumination source to register the phase-shift DHM
    # dxy - Pixel dimensions of the camera sensor used for recording the hologram

    # read dimensions of the holograms
    holo1 = np.array(holo1)
    holo2 = np.array(holo2)
    M, N = holo1.shape

    # The spatial filtering process is executed
    holo_filter1, posi_H1, mask = ut.spatialFilter(holo1, M, N)
    holo_filter2 = ut.spatialFilter2(holo2, M, N, mask)
    if posi_H1[0] > posi_H1[2]:
        posi_H1[0], posi_H1[2] = posi_H1[2], posi_H1[0]
    posi_H2 = posi_H1

    # holo_filter_display = ut.intensity(holo_filter2, True)
    # ut.imageShow(holo_filter_display, 'Filtered Fourier spectrum')

    # create a 3D matrix
    FTHolo = np.zeros((holo1.shape[0], holo2.shape[1], 2), dtype=holo_filter1.dtype)
    FTHolo[:, :, 0] = holo_filter1
    FTHolo[:, :, 1] = holo_filter2

    # apply for loop
    step = 0.025  # rad
    angles = np.arange(0.0, 2 * np.pi + step, step)

    metric_temp1 = 10000
    metric_temp2 = 10000
    array_metric1 = []
    array_metric2 = []
    for cont, phase in enumerate(angles):
        theta_final = [0, phase]
        Gdemod = ut.demComp2SIDHM(theta_final, FTHolo)
        Gplus_demod = Gdemod[:, :, 0]
        Gminus_demod = Gdemod[:, :, 1]

        factor_Gplus = np.abs(Gplus_demod[posi_H1[1], posi_H1[0]]) + np.abs(Gplus_demod[posi_H1[3], posi_H1[2]])
        factor_Gminus = np.abs(Gminus_demod[posi_H2[1], posi_H2[0]]) + np.abs(Gminus_demod[posi_H2[3], posi_H2[2]])

        metric_h1 = np.abs(Gplus_demod[posi_H1[3], posi_H1[2]]) / factor_Gplus
        metric_h2 = np.abs(Gminus_demod[posi_H2[1], posi_H2[0]]) / factor_Gminus
        array_metric1.append(metric_h1)
        array_metric2.append(metric_h2)

        if metric_h1 < metric_temp1:
            Gplus_demod_out = Gplus_demod
            theta_out_Gplus = [0, phase]
            metric_temp1 = metric_h1


        if metric_h2 < metric_temp2:
            Gminus_demod_out = Gminus_demod
            theta_out_Gminus = [0, phase]
            metric_temp2 = metric_h2
    
    #print(theta_out_Gplus, theta_out_Gminus)

    # uncommnet this line if you want to visualize the metris plot profiles
    '''
    plt.plot(angles, array_metric1, label='metric 1')
    plt.plot(angles, array_metric2, label='metric 2')
    plt.xlabel('angles')
    plt.ylabel('values')
    plt.legend()
    plt.show()
    '''

    # Uncomment the lines below to visualize Gplus and Gminus after demodulation
    '''
    Gplus_demod_display = ut.intensity(Gplus_demod_out, True)
    ut.imageShow(Gplus_demod_display, 'Gplus_demod')

    Gminus_demod_display = ut.intensity(Gminus_demod_out, True)
    ut.imageShow(Gminus_demod_display, 'Gminus_demod')
    '''

    # compensation for gplus
    gplus = ut.ift(Gplus_demod_out)
    gplus = gplus * -1
    output_gplus = ut.compensation(gplus, wavelength, dxy, posi_H1[0], posi_H1[1], 2, 0.1)
    # phase = ut.phase(output_gplus)
    # ut.imageShow(phase, 'phase g+')

    # compensation for gminus
    gminus = ut.ift(Gminus_demod_out)
    output_gminus = ut.compensation(gminus, wavelength, dxy, posi_H2[2], posi_H2[3], 2, 0.1)
    #phase = ut.phase(output_gminus)
    #ut.imageShow(phase, 'phase g-')

    # High resolution phase image
    FT_gplus = ut.ft(output_gplus)
    nor_FT_gplus = (FT_gplus - np.min(FT_gplus)) / (np.max(FT_gplus) + np.min(FT_gplus))

    FT_gminus = ut.ft(output_gminus)
    nor_FT_gminus = (FT_gminus - np.min(FT_gminus)) / (np.max(FT_gminus) + np.min(FT_gminus))

    highResolution = nor_FT_gplus + nor_FT_gminus
    #ut.imageShow(ut.intensity(highResolution, False), 'FT High Resolution')

    highResolution = ut.ift(highResolution)

    return highResolution

def siDHM_heuristic(holo1, holo2, wavelength, dxy):
    # Function to retrieve phase maps of off-axis SI-DHM via for loop.
    # inputs:
    # holo1 - [1 of 2] phase-shift holograms
    # holo2 - [2 of 2] phase-shift holograms
    # wavelength - Wavelength of the illumination source to register the phase-shift DHM
    # dxy - Pixel dimensions of the camera sensor used for recording the hologram

    # read dimensions of the holograms
    holo1 = np.array(holo1)
    holo2 = np.array(holo2)
    M, N = holo1.shape

    # The spatial filtering process is executed
    holo_filter1, posi_H1, mask = ut.spatialFilter(holo1, M, N)
    holo_filter2 = ut.spatialFilter2(holo2, M, N, mask)
    if posi_H1[0] > posi_H1[2]:
        posi_H1[0], posi_H1[2] = posi_H1[2], posi_H1[0]
    posi_H2 = posi_H1

    # create a 3D matrix
    FTHolo = np.zeros((holo1.shape[0], holo2.shape[1], 2), dtype=holo_filter1.dtype)
    FTHolo[:, :, 0] = holo_filter1
    FTHolo[:, :, 1] = holo_filter2

    # Define parameters for applying the heuristic search
    step = 0.025
    angles = np.arange(0.0, 2 * np.pi + step, step)
    theta_seeds = [0, random.choice(angles)]

    # Optimization options
    options = {
        'maxiter': 100,  # Maximum number of iterations
        'gtol': 1e-3  # Gradient norm must be less than gtol before successful termination
    }

    # minimization for Gplus
    result = minimize(
        fun=lambda t: ut.costFunction_SIDHM(t, FTHolo, posi_H1),
        x0=theta_seeds, method='BFGS', options=options
    )

    # Extract results
    theta = result.x
    theta_out = [theta[0], theta[1]]
    print(f'angles: {theta_out}')

    Gdemod_out = ut.demComp2SIDHM(theta_out, FTHolo)
    Gplus_demod_out = Gdemod_out[:, :, 0]

    gplus = ut.ift(Gplus_demod_out)
    gplus = gplus * -1
    output_gplus = ut.compensation(gplus, wavelength, dxy, posi_H1[0], posi_H1[1], 2, 0.1)
    # phase = ut.phase(output_gplus)
    # ut.imageShow(phase, 'phase g+')

    # minimization for Gminus
    result = minimize(
        fun=lambda t: ut.costFunction_SIDHM_II(t, FTHolo, posi_H1),
        x0=theta_seeds, method='BFGS', options=options
    )

    # Extract results
    theta = result.x
    theta_out = [theta[0], theta[1]]
    print(f'angles: {theta_out}')

    Gdemod_out = ut.demComp2SIDHM(theta_out, FTHolo)
    Gminus_demod_out = Gdemod_out[:, :, 1]

    gminus = ut.ift(Gminus_demod_out)
    output_gminus = ut.compensation(gminus, wavelength, dxy, posi_H1[2], posi_H1[3], 2, 0.1)
    # phase = ut.phase(output_gminus)
    # ut.imageShow(phase, 'phase g-')

    # High resolution phase image
    FT_gplus = ut.ft(output_gplus)
    nor_FT_gplus = (FT_gplus - np.min(FT_gplus)) / (np.max(FT_gplus) + np.min(FT_gplus))

    FT_gminus = ut.ft(output_gminus)
    nor_FT_gminus = (FT_gminus - np.min(FT_gminus)) / (np.max(FT_gminus) + np.min(FT_gminus))

    highResolution = nor_FT_gplus + nor_FT_gminus
    # ut.imageShow(ut.intensity(highResolution, False), 'FT High Resolution')

    highResolution = ut.ift(highResolution)

    return highResolution
