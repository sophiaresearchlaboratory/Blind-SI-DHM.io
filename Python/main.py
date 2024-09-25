


# import libraries
import utilities as ut
import phaseSI_DHM as SIDHM

# Lines to read the holograms
holo1 = ut.imageRead('C:/Users/racastaneq/PycharmProjects/SI-DHM/samples/Experimental1.bmp')
ut.imageShow(holo1, 'Hologram 1')

holo2 = ut.imageRead('C:/Users/racastaneq/PycharmProjects/SI-DHM/samples/Experimental2.bmp')
ut.imageShow(holo2, 'Hologram 2')

# Lines to load the input parameters (same units)
wavelength = 0.632
dxy = 3.75

# Lines to implement the generalized SI-DHM framework using for loop approach
# output = SIDHM.siDHM_forloop(holo1, holo2, wavelength, dxy)
# phase = ut.phase(output)
# ut.imageShow(phase, 'Phase')


# Lines to implement the generalized SI-DHM framework using heuristic approach
output = SIDHM.siDHM_heuristic(holo1, holo2, wavelength, dxy)
phase = ut.phase(output)
ut.imageShow(phase, 'Phase High resolution')

