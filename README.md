# Blind-SI-DHM Generalized Reconstruction Framework for Structured Illumination in Digital Holographic Microscopy (SI-DHM)

## Project Overview
This project focuses on implementing a generalized reconstruction framework for Structured Illumination Digital Holographic Microscopy (SI-DHM). 
In SI-DHM, high spatial frequencies of microscopic samples are encoded into the +1 diffraction term, which allows for super-resolved phase images through computational reconstruction.
The accuracy of these phase maps depends on the proper computational techniques used, particularly in demodulating the super-resolved components of the +1 diffraction term.

Our framework introduces an automatic demodulation process for the laterally shifted object spectrum, compensating for linear phase terms without requiring prior information about phase shifts. This approach enhances the efficiency and accuracy of super-resolved phase image reconstruction, leading to higher-quality microscopic phase images.

## Features
- Super-resolved Phase Image Reconstruction: Allows the retrieval of high spatial frequency components from the +1 diffraction term.
- Automatic Demodulation: No prior phase shift information is required for demodulation of laterally shifted object spectrums.
- Compensation for Linear Phase Terms: Ensures that linear phase terms are properly accounted for, improving the reconstruction quality.
- Generalized Framework: Flexible and adaptable to different SI-DHM setups.

Requirements
To run this project, you will need the following:

MATLAB R2023a or later (or Python if using an alternative version of the code)

![Logo](/Images/logo_OIRL.png)
