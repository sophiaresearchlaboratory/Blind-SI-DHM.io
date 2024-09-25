
# import lybraries
import numpy as np
import math
from PIL import Image, ImageOps
from matplotlib import pyplot as plt
import sys


def spatialFilter(holo1, M, N):
    # Function to create a spatial filter with a rectangle mask in holograms with structure illumination
    # Inputs:
    # holo1 - hologram to be filtered
    # M,N - hologram dimensions

    print("Spatial filtering process started.....")
    ft_holo = ft(holo1)

    # Finding the max peaks for +1 order in I or II quadrant
    height_half = M // 2
    mask = np.zeros((M, N))
    mask[0:height_half - 1, 0:N] = 1
    ft_holo = ft_holo * mask

    abs_fft = np.abs(ft_holo)
    maximum1 = np.amax(abs_fft)
    fy_Gminus, fx_Gminus = np.where(abs_fft == maximum1)

    ft_holo[fy_Gminus, fx_Gminus] = 0
    maximum2 = np.amax(ft_holo)
    fy_Gplus, fx_Gplus = np.where(ft_holo == maximum2)

    print(f'Pixel values Gplus fx: {fx_Gplus} y fy: {fy_Gplus}')
    print(f'Pixel values Gminus fx: {fx_Gminus} y fy: {fy_Gminus}')

    # create the circular mask
    fx_Gplus = fx_Gplus[0]; fy_Gplus = fy_Gplus[0]; fx_Gminus = fx_Gminus[0]; fy_Gminus = fy_Gminus[0]
    center_x = (fx_Gplus + fx_Gminus) / 2
    center_y = (fy_Gplus + fy_Gminus) / 2
    radius = 10 + np.sqrt((fx_Gplus - fx_Gminus) ** 2 + (fy_Gplus - fy_Gminus) ** 2) / 2
    mask = np.zeros_like(abs_fft)
    y, x = np.ogrid[:mask.shape[0], :mask.shape[1]]
    distance_from_center = np.sqrt((x - center_x) ** 2 + (y - center_y) ** 2)
    mask[distance_from_center <= radius] = 1

    # apply the filter
    holo_filter = ft_holo * mask
    array_fxfy = [fx_Gplus, fy_Gplus, fx_Gminus, fy_Gminus]
    # holo_filter_display = intensity(holo_filter, True)
    # imageShow(holo_filter_display, 'Filter Fourier spectrum')

    print("Spatial filtering process finished.")

    return holo_filter, array_fxfy, mask


def spatialFilter2(holo2, M, N, mask):
    # Function to create a spatial filter with a rectangle mask in holograms with structure illumination
    # Inputs:
    # holo2 - hologram2 to be filtered
    # M,N - hologram dimensions
    ft_holo = ft(holo2)
    holo_filter = ft_holo * mask

    return holo_filter


def imageRead(namefile):
    # Function to read an image file from the disk
    # inputs:
    # namefile - direction image to read
    imagen = Image.open(namefile)
    loadImage = ImageOps.grayscale(imagen)

    return loadImage


def imageShow(image, title):
    # Function to display an Image
    # inputs:
    # image - The image to show
    # title - Title of the displayed image
    plt.imshow(image, cmap='gray'), plt.title(title)
    plt.show()

    return


def amplitude(complexField, log):
    # Function to compute the amplitude of a given complex field
    # inputs:
    # complexField - The input complex field to compute the amplitude
    # log - boolean variable to determine if a log representation is applied
    out = np.abs(complexField)

    if log == True:
        out = 20 * np.log(out)

    return out


def intensity(complexField, log):
    # Function to compute the Intensity of a given complex field
    # inputs:
    # complexField - The input complex field to compute the intensity
    # log - boolean variable to determine if a log representation is applied
    out = np.abs(complexField)
    out = out * out

    if log == True:
        out = 20 * np.log(out)
        out[out == np.inf] = 0
        out[out == -np.inf] = 0

    return out


def phase(complexField):
    # Function to compute the phase of a given complex field
    # inputs:
    # complexField - The input complex field to compute the phase
    out = np.angle(complexField)

    return out


def ft(field):
    # Function to compute the Fourier Transform
    # inputs:
    # field - The input to compute the Fourier Transform
    ft = np.fft.fft2(field)
    ft = np.fft.ifftshift(ft)

    return ft


def ift(field):
    # Function to compute the Inverse Fourier Transform
    # inputs:
    # field - The input to compute the Inverse Fourier Transform
    ift = np.fft.ifftshift(field)
    ift = np.fft.ifft2(ift)
    return ift


def imgInfo(img):
    # Function to get image information
    # inputs:
    # img - The input img to get the information
    width, height = img.size
    print(f"Image size: {width} x {height} pixels")

    return width, height


# Function to save an Image
def saveImg(sample, name):
    # inputs:
    # sample - size image Y
    # name - name image
    image_data = ((sample - np.min(sample)) / (np.max(sample) - np.min(sample)) * 255)
    image = Image.fromarray(image_data.astype(np.uint8))
    image.save(name, format='PNG')

    return


def demComp2SIDHM(theta, H):
    # Function to get demodulation
    # inputs:
    # theta - Array 2D
    # H - Filter Fourier transform of holograms
    X, Y, no = H.shape
    D = np.zeros((X, Y, no), dtype=complex)

    M = 0.5 * np.array([[np.exp(1j * theta[0]), np.exp(-1j * theta[0])],
                        [np.exp(1j * theta[1]), np.exp(-1j * theta[1])]])

    Minv = np.linalg.pinv(M)
    D[:, :, 0] = Minv[0, 0] * H[:, :, 0] + Minv[0, 1] * H[:, :, 1]
    D[:, :, 1] = Minv[1, 0] * H[:, :, 0] + Minv[1, 1] * H[:, :, 1]

    return D


def compensation(field, wavelength, dxy, fx, fy, s, step):
    # Function to compensate phase maps of off-axis DHM via the fast ROI search algorithm.
    # Inputs:
    # field - demodulated G+ or G- to be compensated
    # wavelength - Wavelength of the illumination source to register the DHM hologram
    # dxy - Pixel dimensions of the camera sensor used for recording the hologram
    # s = 5 and step = 0.2

    print("Phase compensation started....")

    # Creating a mesh_grid to operate in world-coordinates
    M, N = field.shape
    x = np.arange(0, N, 1)
    y = np.arange(0, M, 1)
    X, Y = np.meshgrid(x - (N / 2), y - (M / 2), indexing='xy')
    fx_0, fy_0 = N / 2, M / 2
    k = (2 * math.pi) / wavelength

    fin = 0
    G_temp = s
    while fin == 0:
        sum_max = 0  # small number for the metric (thresholding)
        arrayY = np.linspace(int(10 * (fy - step * G_temp)), int(10 * (fy + step * G_temp)), int(10 * step))
        arrayX = np.linspace(int(10 * (fx - step * G_temp)), int(10 * (fx + step * G_temp)), int(10 * step))

        for fx_temp in arrayX:
            for fy_temp in arrayY:
                fx_tmp, fy_tmp = fx_temp / 10, fy_temp / 10

                # Digital reference wave
                theta_x = math.asin((fx_0 - fx_tmp) * wavelength / (N * dxy))
                theta_y = math.asin((fy_0 - fy_tmp) * wavelength / (M * dxy))
                ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * dxy) + (math.sin(theta_y) * Y * dxy)))

                # Compensation of the tilting angle for the off-axis acquisition
                reconstruction = field * ref_wave
                phase = np.angle(reconstruction)

                # Thresholding process
                minVal = np.amin(phase)
                maxVal = np.amax(phase)
                phase_sca = (phase - minVal) / (maxVal - minVal)
                binary_phase = (phase_sca > 0.2)

                # Applying the summation and thresholding metric
                sum = np.sum(np.sum(binary_phase))
                if (sum > sum_max):
                    x_max_out = fx_tmp
                    y_max_out = fy_tmp
                    sum_max = sum
                # print (sum)

        G_temp = G_temp - 1

        if x_max_out == fx and y_max_out == fy:
            fin = 1;

        fx = x_max_out
        fy = y_max_out

    # Retrieving the best reconstruction (compensated phase)
    theta_x = math.asin((fx_0 - x_max_out) * wavelength / (N * dxy))
    theta_y = math.asin((fy_0 - y_max_out) * wavelength / (M * dxy))
    ref_wave = np.exp(1j * k * ((math.sin(theta_x) * X * dxy) + (math.sin(theta_y) * Y * dxy)))
    comp_phase = field * ref_wave

    print("Phase compensation finished.")

    return comp_phase


def costFunction_SIDHM(theta, FTHolo, posi_H1):
    # Function to compute the demodulation in the Fourier domain
    # Inputs:
    # theta - demodulated G+ or G- to be compensated
    # wavelength - Wavelength of the illumination source to register the DHM hologram
    # dxy - Pixel dimensions of the camera sensor used for recording the hologram
    # s = 5 and step = 0.2
    Gdemod = demComp2SIDHM(theta, FTHolo)

    # Extract the first layer in the third dimension
    Gplus_demod = Gdemod[:, :, 0]

    # Calculate maxL and maxR using provided indices
    factor_Gplus = np.abs(Gplus_demod[posi_H1[1], posi_H1[0]]) + np.abs(Gplus_demod[posi_H1[3], posi_H1[2]])
    J = np.abs(Gplus_demod[posi_H1[3], posi_H1[2]]) / factor_Gplus

    return J


def costFunction_SIDHM_II(theta, FTHolo, posi_H1):
    # Function to compute the demodulation in the Fourier domain
    # Inputs:
    # theta - demodulated G+ or G- to be compensated
    # wavelength - Wavelength of the illumination source to register the DHM hologram
    # dxy - Pixel dimensions of the camera sensor used for recording the hologram
    # s = 5 and step = 0.2
    Gdemod = demComp2SIDHM(theta, FTHolo)

    # Extract the first layer in the third dimension
    Gplus_demod = Gdemod[:, :, 1]

    # Calculate maxL and maxR using provided indices
    factor_Gplus = np.abs(Gplus_demod[posi_H1[1], posi_H1[0]]) + np.abs(Gplus_demod[posi_H1[3], posi_H1[2]])
    J = np.abs(Gplus_demod[posi_H1[1], posi_H1[0]]) / factor_Gplus

    return J