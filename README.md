# Generalized computational framework for Structured Illumination Digital Holography
HOLAA ANA Y SOPHIA
The use of a structured illumination (SI) pattern to illuminate a target complex object in digital holographic microscopy (DHM) allows for encoding high spatial frequencies of the target microscopic sample into the +1 diffraction term, reconstructing super-resolved phase images after the proper computational reconstruction. The accuracy of these phase maps depends on the correct demodulation of the high-resolution frequencies within the +1 term and the proper compensation of the interference angle between the object and reference waves of the DHM system. This work presents a generalized reconstruction framework for SI-DHM that automatically demodulates the two laterally-shifted object spectrums and compensates for the linear phase term of the optical interference with minimum input from the user (i.e., two recorded holograms with a phase shift of the SI pattern, the source’s wavelength and the pixel size of the sensor). 

## Generalized SI-DHM framework 
The block diagram of the generalized SI-DHM framework contains six stages. The first stage is focused on define the input parameters, which are the two recorded holograms (h<sub>i</sub>), the source wavelength (λ), and the pixel size of the sensor (Δxy). The second one is focused on estimating the Fourier spectra of the recorded holograms (H<sub>i</sub>). The third one corresponds to the spatial filtering of the hologram spectrum to select the +1 term (H+1,i). The fourth and fifth steps are the demodulation of the laterally-shifted object spectrum encoded in the +1 term (i.e., G<sub>+</sub> and G<sub>-</sub>) and its proper centering in the frequency domain, respectively. Finally, the last step is the combination of both centered object spectra to generate the super-resolved phase image.


![Flowchart](/Images/flow.png)

## Methods of use
We have developed two implementations of the proposed generalized SI-DHM framework. One version minimizes the cost functions using iterative loops, while the other uses a heuristic approach. Both implementations are available for download in MATLAB and Python


<div class="table_component" role="region" tabindex="0" style="width:100%">
<table style="border: 1rem;margin-left: auto; margin-right: auto;">
    <thead>
        <tr>
            <th></th>
            <th>Iterative loop</th>
            <th>Heuristic approach</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td><b>Matlab</b></td>
            <td><a href="https://github.com/sophiaresearchlaboratory/Blind-SI-DHM.io/tree/main/Matlab" download>Click to download</a></td>
            <td><a href="https://github.com/sophiaresearchlaboratory/Blind-SI-DHM.io/tree/main/Matlab/heuristic" download>Click to download</a></td>
        </tr>
        <tr>
            <td><b>Python</b></td>
            <td><a href="https://drive.google.com/drive/folders/18gG1earq-w9dI_Nl5ZCrUP7y6cUk4riR?usp=drive_link" download>Click to download</a></td>
            <td><a href="https://drive.google.com/drive/folders/18gG1earq-w9dI_Nl5ZCrUP7y6cUk4riR?usp=drive_link" download>Click to download</a></td>
        </tr>
    </tbody>
</table>
</div>


## Download Samples
<a href="https://github.com/sophiaresearchlaboratory/Blind-SI-DHM.io/tree/main/Samples/start">Test Start</a> 
Parameter: wavelength 632.8 nm and pixel size 3.75 um

<a href="https://github.com/sophiaresearchlaboratory/Blind-SI-DHM.io/tree/main/Samples/USAF">Test USAF</a> 
Parameters: wavelength 532 nm and pixel size 3.75 um

## Example of Matlab
| ![Descripción de la imagen](Images/Matlab.png) |

## Example of Python
| ![Descripción de la imagen](Images/python.png) |

## Funding
This project was funded by Vicerrectoría de Ciencia, Tecnología e Innovacion, the Fundamental Sciences area at Universidad EAFIT, and the National Science Foundation (NSF) (grant number 2042563 and 2404769).

## Credits
-	MATLAB (2020). version 7.10.0 (R2020a). Natick, Massachusetts: The MathWorks Inc.
*	Python Software Foundation. 

## Citation
If using the generalized SI-DHM framework for publication, please kindly cite the following: S. Obando-Vasquez, R. Castaneda, R. Restrepo, C. Trujilo, and A. Doblas, **“Generalized computational framework for phase image reconstruction in structured illumination digital holographic microscopy,”** Optics Express, under review (2024). 

## Support or Contact
<table>
    <tr>
        <td><b>Researcher</b></td>
        <td><b>e-mail</b></td>
        <td><b>Google Shocolar</b></td>
    </tr>
    <tr>
        <td>Sofia Obando</td>
        <td>sobandovasquez@umassd.edu</td>
        <td><a href="https://scholar.google.com/citations?user=01-dFrAAAAAJ&hl=es">Google Scholar</a></td>
    </tr>
    <tr>
        <td>Raul Castaneda</td>
        <td>racastaneq@eafit.edu.co</td>
        <td><a href="https://scholar.google.com.co/citations?user=RBtkL1oAAAAJ&hl=en">Google Scholar</a></td>
    </tr>
    <tr>
        <td>Rene Restrepo</td>
        <td>rrestre6@eafit.edu.co</td>
        <td><a href="https://scholar.google.com.co/citations?hl=en&user=nLcBLUwAAAAJ&view_op=list_works&sortby=pubdate">Google Scholar</a></td>
    </tr>
    <tr>
        <td>Carlos Trujillo*</td>
        <td>catrujilla@eafit.edu.co</td>
        <td><a href="https://scholar.google.com.co/citations?user=BKVrl2gAAAAJ&hl=en">Google Scholar</a></td>
    </tr>
    <tr>
        <td>Ana Doblas*</td>
        <td>adoblas@umassd.edu</td>
        <td><a href="https://scholar.google.com.co/citations?user=PvvDEMYAAAAJ&hl=en">Google Scholar</a></td>
    </tr>
</table>
The Principal Investigators (PIs) of the project are Drs. Carlos Trujillo, and Ana Doblas


<img src="https://raw.githubusercontent.com/sophiaresearchlaboratory/Blind-SI-DHM.io/refs/heads/main/Images/logo_OIRL.png" width="400" height="87" align="left"> <img src="https://raw.githubusercontent.com/sophiaresearchlaboratory/Blind-SI-DHM.io/refs/heads/main/Images/SOPHIA_RG.png" width="426" height="87" align="right">
